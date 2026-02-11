from __future__ import annotations
from pathlib import Path
from typing import Optional, Tuple
import numpy as np
from netCDF4 import Dataset

from xy2ll import xy2ll
from ll2xy import ll2xy
try:
    from InterpFromGridToMesh import InterpFromGridToMesh  # returns (Fq, in_mesh)
except ImportError:
    InterpFromGridToMesh = None

# interp_searise.py
# robust_interp_searise.py
import socket
import numpy as np



def _default_searise_path(hemisphere: int) -> str:
    host = socket.gethostname()
    if host in {"ronne"}:
        return "/home/ModelData/SeaRISE/Greenland_5km_dev1.2.nc" if hemisphere==+1 else "/home/ModelData/SeaRISE/Antarctica_5km_dev1.0.nc"
    if host in {"thwaites","larsen","murdo","astrid"}:
        return "/u/astrid-r1b/ModelData/SeaRISE/Greenland5km_v1.2/Greenland_5km_dev1.2.nc" if hemisphere==+1 else "/u/astrid-r1b/ModelData/SeaRISE/Antarctica5km_shelves_v1.0/Antarctica_5km_dev1.0.nc"
    if host in {"totten"}:
        return "/totten_1/ModelData/SeaRISE/Greenland_5km_dev1.2.nc" if hemisphere==+1 else "/totten_1/ModelData/SeaRISE/Antarctica_5km_dev1.0.nc"
    if host in {"nri-085597"}:
        return "aq1_01_20.nc" if hemisphere==+1 else "aq1_01_20.nc"
    raise RuntimeError(f"hostname '{host}' not supported; pass ncfile=")


def _prepare_grid(xdata, ydata, raw):
    """Return (x, y, data_yx) so that:
       - x.shape=(nx,), y.shape=(ny,)
       - data_yx.shape=(ny, nx) matches InterpFromGrid(y,x) convention
       - axes strictly increasing; data reindexed if we reversed an axis
    """
    x = np.asarray(xdata).astype(float).ravel()
    y = np.asarray(ydata).astype(float).ravel()
    A = np.asarray(raw).astype(float)

    # Handle potential trailing dims (e.g., (ny, nx, 1))
    while A.ndim > 2 and 1 in A.shape:
        A = np.squeeze(A)

    # Decide whether to transpose based on shapes
    if A.shape == (y.size, x.size):
        data_yx = A
    elif A.shape == (x.size, y.size):
        data_yx = A.T
    else:
        # Fall back: if one dim matches x and the other matches y, permute
        if x.size in A.shape and y.size in A.shape and A.ndim == 2:
            ix = A.shape.index(y.size)
            jx = A.shape.index(x.size)
            data_yx = A if (ix == 0 and jx == 1) else A.transpose()
        else:
            raise ValueError(f"Unexpected data shape {A.shape}; x={x.size}, y={y.size}")

    # Ensure axes are increasing, and permute data accordingly
    if np.any(np.diff(x) < 0):
        x = x[::-1]
        data_yx = data_yx[:, ::-1]
    if np.any(np.diff(y) < 0):
        y = y[::-1]
        data_yx = data_yx[::-1, :]

    # Heuristic: SeaRISE coords are meters; if values look like km, scale up
    # (e.g., max ~3500 suggests km for Antarctica domain)
    if np.nanmax(np.abs(x)) < 10000 and np.nanmax(np.abs(y)) < 10000:
        x = x * 1000.0
        y = y * 1000.0

    return x, y, data_yx


def interp_searise(X, Y, varname, hemisphere=+1, ncfile=None, verbose=False):
    if ncfile is None:
        ncfile = _default_searise_path(hemisphere)

    X = np.asarray(X)
    Y = np.asarray(Y)
    xq = X.ravel().astype(float)
    yq = Y.ravel().astype(float)

    # Projection logic: matches your MATLAB exactly
    if hemisphere == +1:
        LAT, LON = xy2ll(xq, yq, +1, 45, 70)
        xproj, yproj = ll2xy(LAT, LON, +1, 39, 71)
    elif hemisphere == -1:
        xproj, yproj = xq, yq
    else:
        raise ValueError("hemisphere must be +1 (Greenland) or -1 (Antarctica)")

    with Dataset(ncfile, "r") as ds:
        xdata = ds.variables["rlon"][:]
        ydata = ds.variables["rlat"][:]
        raw = ds.variables[varname][:]

        A = np.asarray(raw)

        # --- NEW: collapse trivial dimensions so RACMO shapes like
        # (nt, 1, ny, nx) become (nt, ny, nx) -----------------------------
        # First, special-case 4D with a singleton 2nd dim (time,1,y,x)
        if A.ndim == 4 and A.shape[1] == 1:
            A = A[:, 0, :, :]  # -> (nt, ny, nx)

        # Then, for anything >=3D, squeeze any remaining size-1 dims
        # but keep time axis intact in typical (nt, ny, nx) case.
        if A.ndim > 2 and 1 in A.shape:
            A = np.squeeze(A)
        # ------------------------------------------------------------------

        # Now interpret shape
        if A.ndim == 3:
            # Assume (nt, ny, nx): take climatological mean over time
            if verbose:
                print(
                    f"interp_searise: '{varname}' is 3D {A.shape}, "
                    "taking nanmean over time axis 0 → 2D field"
                )
            raw2d = np.nanmean(A, axis=0)
        elif A.ndim == 2:
            raw2d = A
        else:
            raise ValueError(
                f"interp_searise: cannot handle '{varname}' with ndim={A.ndim}, "
                f"shape={A.shape}"
            )

    # Auto-fix orientation & axis order for (x, y, data_yx)
    xg, yg, data_yx = _prepare_grid(xdata, ydata, raw2d)

    # Interp (nearest for LandMask, linear otherwise) – identical to MATLAB call
    if varname.lower() == "landmask":
        vals = InterpFromGridToMesh(xg, yg, data_yx, xproj, yproj, 0, "nearest")
    else:
        vals = InterpFromGridToMesh(xg, yg, data_yx, xproj, yproj, 0)

    return vals.reshape(X.shape)

