from __future__ import annotations
from pathlib import Path
import numpy as np
import torch
import torch.nn as nn

DEFAULT_WEIGHTS_PATH = str(Path(__file__).resolve().parent / "friction_emulator.pt")

X_FLOOR = np.array([1.0e-30, 1.0e-12], dtype=np.float32)
Y_FLOOR = np.array([1.0e-30], dtype=np.float32)
X_MEAN = np.array([15.853013060672888, -14.369379188536289], dtype=np.float32)
X_STD = np.array([1.2777508831457713, 1.0928669173617296], dtype=np.float32)
Y_MEAN = np.array([25.432599186435276], dtype=np.float32)
Y_STD = np.array([1.4568046734404581], dtype=np.float32)
NORMALIZATION = "log"
VMAG_MIN = 5.0e-8

_MODEL = None
_DEVICE = None
_NORMALIZATION = NORMALIZATION
_VMAG_MIN = VMAG_MIN
_X_FLOOR_T = None
_X_MEAN_T = None
_X_STD_T = None
_Y_MEAN_T = None
_Y_STD_T = None

class FrictionMLP(nn.Module):# {{{
    def __init__(self, in_dim: int = 2, h1: int = 64, h2: int = 64, out_dim: int = 1):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, h1),
            nn.ReLU(),
            nn.Linear(h1, h2),
            nn.ReLU(),
            nn.Linear(h2, out_dim),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)# }}}
def init_model(weights_path: str = DEFAULT_WEIGHTS_PATH, device: str = "auto") -> None:# {{{
    global _MODEL, _DEVICE, _NORMALIZATION, _VMAG_MIN, _X_FLOOR_T, _X_MEAN_T, _X_STD_T, _Y_MEAN_T, _Y_STD_T

    ckpt_path = Path(weights_path)
    if not ckpt_path.exists():
        raise FileNotFoundError(f"Friction emulator checkpoint not found: {ckpt_path}")

    checkpoint = torch.load(str(ckpt_path), map_location="cpu", weights_only=False)
    in_dim = int(checkpoint["in_dim"])
    h1 = int(checkpoint["h1"])
    h2 = int(checkpoint["h2"])
    out_dim = int(checkpoint["out_dim"])

    if device == "auto":
        resolved_device = "cuda" if torch.cuda.is_available() else "cpu"
    else:
        resolved_device = device

    model = FrictionMLP(in_dim=in_dim, h1=h1, h2=h2, out_dim=out_dim)
    model.load_state_dict(checkpoint["state_dict"])
    model.to(resolved_device)
    model.eval()

    normalization = str(checkpoint.get("transform", checkpoint.get("normalization", NORMALIZATION)))
    if normalization == "none":
        normalization = "raw"
    x_floor = np.asarray(checkpoint.get("x_floor", X_FLOOR), dtype=np.float32)
    x_mean = np.asarray(checkpoint.get("x_mean", X_MEAN), dtype=np.float32)
    x_std = np.asarray(checkpoint.get("x_std", X_STD), dtype=np.float32)
    y_mean = np.asarray(checkpoint.get("y_mean", Y_MEAN), dtype=np.float32)
    y_std = np.asarray(checkpoint.get("y_std", Y_STD), dtype=np.float32)
    vmag_min = float(checkpoint.get("min_vmag", VMAG_MIN))

    _MODEL = model
    _DEVICE = resolved_device
    _NORMALIZATION = normalization
    _VMAG_MIN = vmag_min
    _X_FLOOR_T = torch.as_tensor(x_floor, dtype=torch.float32, device=resolved_device)
    _X_MEAN_T = torch.as_tensor(x_mean, dtype=torch.float32, device=resolved_device)
    _X_STD_T = torch.as_tensor(x_std, dtype=torch.float32, device=resolved_device)
    _Y_MEAN_T = torch.as_tensor(y_mean, dtype=torch.float32, device=resolved_device)
    _Y_STD_T = torch.as_tensor(y_std, dtype=torch.float32, device=resolved_device)
    print(
        f"Friction emulator initialized on device: {resolved_device}, "
        f"normalization: {normalization}, min_vmag: {vmag_min:.6e}"
    )# }}}


def _transform_features(feats_t: torch.Tensor) -> torch.Tensor:
    if _NORMALIZATION in ("raw", "standard"):
        return feats_t
    if _NORMALIZATION == "sqrt":
        return torch.sqrt(torch.clamp_min(feats_t, 0.0))
    if _NORMALIZATION == "mixed":
        c2 = feats_t[:, 0:1]
        vmag = torch.log(torch.maximum(feats_t[:, 1:2], _X_FLOOR_T[1]))
        return torch.cat((c2, vmag), dim=1)
    if _NORMALIZATION == "log":
        return torch.log(torch.maximum(feats_t, _X_FLOOR_T))
    raise RuntimeError(f"Unsupported normalization mode: {_NORMALIZATION}")


def _inverse_transform_alpha2(pred_trans: torch.Tensor) -> torch.Tensor:
    if _NORMALIZATION in ("raw", "standard"):
        return torch.clamp_min(pred_trans, 0.0)
    if _NORMALIZATION == "sqrt":
        return torch.square(torch.clamp_min(pred_trans, 0.0))
    if _NORMALIZATION in ("log", "mixed"):
        return torch.exp(pred_trans)
    raise RuntimeError(f"Unsupported normalization mode: {_NORMALIZATION}")


@torch.no_grad()
def predict_alpha2_np(feats, *, dtype="float64"):# {{{
    if _MODEL is None:
        raise RuntimeError("Friction emulator is not initialized")

    feats_np = np.asarray(feats, dtype=np.float32)
    if feats_np.ndim == 1:
        feats_np = feats_np.reshape(1, -1)
    if feats_np.shape[1] != 2:
        raise ValueError(f"Expected input shape (*, 2), got {feats_np.shape}")

    feats_np = feats_np.copy()
    feats_np[:, 1] = np.maximum(feats_np[:, 1], np.float32(_VMAG_MIN))

    feats_t = torch.as_tensor(feats_np, dtype=torch.float32, device=_DEVICE)
    feats_trans = _transform_features(feats_t)
    feats_norm = (feats_trans - _X_MEAN_T) / _X_STD_T
    pred_norm = _MODEL(feats_norm)
    pred_trans = pred_norm * _Y_STD_T + _Y_MEAN_T
    pred_raw = _inverse_transform_alpha2(pred_trans)
    pred_raw = pred_raw.detach().cpu().contiguous()
    pred_raw = pred_raw.to(getattr(torch, dtype)) if isinstance(dtype, str) else pred_raw.to(dtype)
    return pred_raw.numpy().copy()# }}}
