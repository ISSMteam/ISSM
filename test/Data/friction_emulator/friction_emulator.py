from __future__ import annotations
from pathlib import Path
import numpy as np
import torch
import torch.nn as nn

DEFAULT_WEIGHTS_PATH = (
    "/home1/10783/yinmin/work/Applications/ISSM/src/c/modules/"
    "FrictionEmulator/trained_models/friction_emulator.pt"
)

FIXED_X_MEAN = np.array([9.05e6, 2.08e-5], dtype=np.float32)
FIXED_X_STD = np.array([6.61e6, 4.67e-5], dtype=np.float32)
FIXED_Y_MEAN = np.array([2.09e11], dtype=np.float32)
FIXED_Y_STD = np.array([1.18e12], dtype=np.float32)

_MODEL = None
_DEVICE = None
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
    global _MODEL, _DEVICE, _X_MEAN_T, _X_STD_T, _Y_MEAN_T, _Y_STD_T

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

    _MODEL = model
    _DEVICE = resolved_device
    _X_MEAN_T = torch.as_tensor(FIXED_X_MEAN, dtype=torch.float32, device=resolved_device)
    _X_STD_T = torch.as_tensor(FIXED_X_STD, dtype=torch.float32, device=resolved_device)
    _Y_MEAN_T = torch.as_tensor(FIXED_Y_MEAN, dtype=torch.float32, device=resolved_device)
    _Y_STD_T = torch.as_tensor(FIXED_Y_STD, dtype=torch.float32, device=resolved_device)
    print(f"Friction emulator initialized on device: {resolved_device}")# }}}
@torch.no_grad()
def predict_alpha2_np(feats, *, dtype="float64"):# {{{
    if _MODEL is None:
        raise RuntimeError("Friction emulator is not initialized")

    feats_np = np.asarray(feats, dtype=np.float32)
    if feats_np.ndim == 1:
        feats_np = feats_np.reshape(1, -1)
    if feats_np.shape[1] != 2:
        raise ValueError(f"Expected input shape (*, 2), got {feats_np.shape}")

    feats_t = torch.as_tensor(feats_np, dtype=torch.float32, device=_DEVICE)
    feats_norm = (feats_t - _X_MEAN_T) / _X_STD_T
    pred_norm = _MODEL(feats_norm)
    pred_raw = pred_norm * _Y_STD_T + _Y_MEAN_T
    pred_raw = pred_raw.detach().cpu().contiguous()
    pred_raw = pred_raw.to(getattr(torch, dtype)) if isinstance(dtype, str) else pred_raw.to(dtype)
    return pred_raw.numpy().copy()# }}}
