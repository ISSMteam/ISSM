"""Daily ERA5 to ISSM SMB inference for the six-predictor Greenland U-Net."""
from datetime import date, datetime, timedelta
from functools import lru_cache
import math
import os
from pathlib import Path

import h5py

import numpy as np
import torch
from torch import nn
from torch.nn import functional as F


def block(cin, cout):
    return nn.Sequential(nn.Conv2d(cin, cout, 3, padding=1), nn.GroupNorm(1, cout), nn.ReLU(),
                         nn.Conv2d(cout, cout, 3, padding=1), nn.GroupNorm(1, cout), nn.ReLU())


class UNet(nn.Module):
    def __init__(self, in_channels, out_channels, base_channels=32):
        super().__init__()
        b = base_channels
        self.enc1, self.enc2 = block(in_channels, b), block(b, 2*b)
        self.bridge = block(2*b, 4*b)
        self.dec2, self.dec1 = block(6*b, 2*b), block(3*b, b)
        self.head = nn.Conv2d(b, out_channels, 1)

    def forward(self, x):
        a = self.enc1(x)
        b = self.enc2(F.max_pool2d(a, 2))
        c = self.bridge(F.max_pool2d(b, 2))
        d = self.dec2(torch.cat([F.interpolate(c, size=b.shape[-2:], mode='bilinear', align_corners=False), b], 1))
        e = self.dec1(torch.cat([F.interpolate(d, size=a.shape[-2:], mode='bilinear', align_corners=False), a], 1))
        return self.head(e)


SECONDS_PER_DAY = 86400
FEATURES = ('t2m', 'mtpr', 'msr', 'msdwlwrf', 'msdwswrf', 'sfcwind', 'heightsel')
TARGETS = ('sf', 'rf', 'ru', 'su')
ERA_UNITS = {
    't2m': 'K', 'mtpr': 'kg m**-2 s**-1', 'msr': 'kg m**-2 s**-1',
    'msdwlwrf': 'W m**-2', 'msdwswrf': 'W m**-2',
    'u10': 'm s**-1', 'v10': 'm s**-1',
}
_state = None


def init_model(checkpoint_path, device='auto'):
    """Load a six-predictor Greenland checkpoint and the adjacent grid and mask."""
    global _state
    path = Path(checkpoint_path).resolve()
    forcing_dir = Path(os.environ['ISSM_UNET_ERA_DIR']).expanduser().resolve()
    start_date = date.fromisoformat(os.environ['ISSM_UNET_FORCING_START_DATE'])
    if not forcing_dir.is_dir():
        raise ValueError(f'ERA forcing directory does not exist: {forcing_dir}')
    if device == 'auto':
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
    device = torch.device(device)
    checkpoint = torch.load(path, map_location='cpu', weights_only=True)
    config = checkpoint['config']
    if tuple(config['features']) != FEATURES or tuple(config['targets']) != TARGETS:
        raise ValueError('Expected the six ERA predictors plus elevation, with SF/RF/RU/SU targets')
    ny, nx = shape = tuple(config['shape'])
    with np.load(path.with_name('region_masks.npz'), allow_pickle=False) as data:
        mask = data['Greenland'].astype(bool)
        reference = data['SH'].astype(np.float32)
        lat, lon = data['LAT'], data['LON']
    with np.load(path.with_name('native_grid.npz'), allow_pickle=False) as data:
        x = data['X14_311'].astype(np.float64) * 1000
        y = data['Y21_560'].astype(np.float64) * 1000
        if not np.array_equal(lat, data['LAT']) or not np.array_equal(lon, data['LON']):
            raise ValueError('Native grid coordinates differ from the training mask grid')
    if (mask.shape != shape or reference.shape != shape or
            not mask.any() or not np.isfinite(reference[mask]).all()):
        raise ValueError('Invalid training mask or reference elevation')
    for axis, size in ((x, nx), (y, ny)):
        if axis.shape != (size,) or not np.isfinite(axis).all() or not (np.diff(axis) > 0).all():
            raise ValueError('Invalid native grid axis')
    stats = {}
    for key, size in (('x_mean', len(FEATURES)), ('x_std', len(FEATURES)),
                      ('y_mean', len(TARGETS)), ('y_std', len(TARGETS))):
        values = np.asarray(config['normalization'][key], dtype=np.float32)
        if (values.shape != (size,) or not np.isfinite(values).all() or
                (key.endswith('_std') and not (values > 0).all())):
            raise ValueError(f'Invalid checkpoint normalization: {key}')
        stats[key] = values[:, None, None]
    args = config['model_args']
    if args['in_channels'] != len(FEATURES) + 1 or args['out_channels'] != len(TARGETS):
        raise ValueError('Checkpoint channel counts do not match its inputs and targets')
    model = UNet(**args)
    model.load_state_dict(checkpoint['state_dict'], strict=True)
    model.to(device).eval()
    _state = dict(model=model, device=device, shape=shape, stats=stats, mask=mask,
                  reference=reference, x=x, y=y, forcing_dir=forcing_dir, start_date=start_date)
    _read_forcing.cache_clear()
    return dict(shape=shape, epoch=checkpoint['epoch'])


def _text(value):
    return value.decode() if isinstance(value, bytes) else value


@lru_cache(maxsize=1)
def _read_forcing(forcing_date):
    """Read and cache a daily remapped ERA plane in the training feature units."""
    state = _state
    path = state['forcing_dir'] / f'ERA5-{forcing_date.year}.nc'
    with h5py.File(path, 'r') as era:
        for name in ('x', 'y'):
            axis = era[name]
            if (_text(axis.attrs['units']) != 'km' or
                    not np.array_equal(axis[:].astype(np.float64) * 1000, state[name])):
                raise ValueError(f'{path}: {name} does not match the training grid')
        time = era['time']
        if _text(time.attrs.get('calendar', 'standard')) not in ('standard', 'gregorian', 'proleptic_gregorian'):
            raise ValueError(f'{path}: unsupported ERA calendar')
        unit, sep, origin = _text(time.attrs['units']).partition(' since ')
        if not sep or unit.lower() not in ('day', 'days'):
            raise ValueError(f'{path}: expected ERA time in days since an epoch')
        try:
            epoch = datetime.fromisoformat(origin)
        except ValueError:
            epoch = datetime.strptime(origin, '%d-%b-%Y %H:%M:%S')
        matches = [i for i, offset in enumerate(time[:])
                   if (epoch + timedelta(days=float(offset))).date() == forcing_date]
        if len(matches) != 1:
            raise ValueError(f'{path}: expected exactly one forcing record for {forcing_date}')

        fields = {}
        for name, units in ERA_UNITS.items():
            field = era[name]
            scales = tuple(tuple(scale.name for scale in dim.values()) for dim in field.dims)
            if field.shape != (len(time), *state['shape']) or scales != (('/time',), ('/y',), ('/x',)):
                raise ValueError(f'{path}: expected {name} on the (time,y,x) grid')
            if _text(field.attrs['units']) != units:
                raise ValueError(f'{path}: expected {name} units {units}')
            raw = field[matches[0]]
            invalid = ~np.isfinite(raw)
            for attr in ('_FillValue', 'missing_value'):
                for marker in np.asarray(field.attrs.get(attr, [])).ravel():
                    invalid |= raw == marker
            values = raw.astype(np.float32) * float(np.asarray(field.attrs.get('scale_factor', 1)).item())
            values += float(np.asarray(field.attrs.get('add_offset', 0)).item())
            values[invalid] = np.nan
            fields[name] = values
    fields['t2m'] -= np.float32(273.15)
    fields['sfcwind'] = np.hypot(fields.pop('u10'), fields.pop('v10'))
    return fields


def predict_smb_np(surface_grid_m, forcing_date, rho_ice, x_grid_m, y_grid_m):
    """Return reconstructed SMB in m ice/s, with reference SH beyond the ISSM mesh."""
    state = _state
    if state is None:
        raise RuntimeError('Call init_model first')
    if isinstance(forcing_date, str):
        forcing_date = date.fromisoformat(forcing_date)
    if not np.isfinite(rho_ice) or rho_ice <= 0:
        raise ValueError('Ice density must be positive and finite')
    if not np.array_equal(x_grid_m, state['x']) or not np.array_equal(y_grid_m, state['y']):
        raise ValueError('ISSM interpolation axes differ from the training grid in metres')
    surface = np.asarray(surface_grid_m, dtype=np.float32)
    if surface.shape != state['shape'] or np.isinf(surface).any():
        raise ValueError('Invalid ISSM surface grid shape or values')
    # C++ uses NaN for uncovered cells. A finite zero is a valid elevation.
    elevation = np.where(np.isnan(surface), state['reference'], surface)
    fields = _read_forcing(forcing_date)
    raw = np.stack([elevation if name == 'heightsel' else fields[name] for name in FEATURES])
    mask = state['mask'] & np.isfinite(raw).all(axis=0)
    if not mask.any():
        raise ValueError(f'No valid U-Net inputs for {forcing_date}')
    stats = state['stats']
    normalized = np.where(mask, (raw - stats['x_mean']) / stats['x_std'], 0)
    model_input = np.concatenate([normalized, mask[None].astype(np.float32)])
    with torch.inference_mode():
        output = state['model'](torch.from_numpy(model_input[None]).to(state['device']))
        components = output[0].cpu().numpy() * stats['y_std'] + stats['y_mean']
    if not np.isfinite(components[:, mask]).all():
        raise ValueError(f'Nonfinite U-Net predictions for {forcing_date}')
    sf, rf, ru, su = components.astype(np.float64)
    smb = (sf + rf - ru - su) / (rho_ice * SECONDS_PER_DAY)
    return np.ascontiguousarray(np.where(mask, smb, 0.0))


def predict_smb_at_time_np(surface_grid_m, elapsed_seconds, rho_ice, x_grid_m, y_grid_m):
    """Map current ISSM elapsed seconds to a Gregorian forcing date."""
    if _state is None:
        raise RuntimeError('Call init_model first')
    if not np.isfinite(elapsed_seconds) or elapsed_seconds < 0:
        raise ValueError('ISSM elapsed time must be finite and nonnegative')
    day = math.floor(elapsed_seconds / SECONDS_PER_DAY + 1e-9)
    forcing_date = _state['start_date'] + timedelta(days=day)
    return predict_smb_np(surface_grid_m, forcing_date, rho_ice, x_grid_m, y_grid_m)
