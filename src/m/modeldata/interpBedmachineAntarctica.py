# interp_bedmachine_antarctica.py
from __future__ import annotations

from pathlib import Path
import os
import numpy as np
from netCDF4 import Dataset

# SciPy is optional but recommended for 'cubic'
try:
    from scipy.interpolate import RegularGridInterpolator, RectBivariateSpline
    _HAVE_SCIPY = True
except Exception:
    _HAVE_SCIPY = False


def interp_bedmachine_antarctica(X, Y, string="bed", method=None, ncdate="v3"):
    """
    Python port of interpBedmachineAntarctica(X,Y,string,method,ncdate)

    Parameters
    ----------
    X, Y : array-like
        Target coordinates (same shape) in the BedMachine (EPSG:3031) grid units (meters).
    string : {"bed","surface","thickness","mask","source",...,"icemask"}
        Variable name in the NetCDF, plus special case "icemask".
    method : {"nearest","linear","cubic"} or None
        Interpolation method. If None: 'nearest' for mask/source, else 'cubic'.
    ncdate : str
        Either a BedMachine version tag (e.g. "v3.5") to resolve common paths,
        or a full path to a specific .nc file.

    Returns
    -------
    output : ndarray
        Interpolated values with the same shape as X and Y.
    """
    if method is None:
        if string in ("mask", "source"):
            method = "nearest"
        else:
            method = "cubic"  # default, like MATLAB

    # Resolve dataset path
    basename = "BedMachineAntarctica"
    if (str(ncdate).endswith(".nc") and Path(ncdate).exists()) or Path(ncdate).is_file():
        ncfile = str(ncdate)
        version_str = Path(ncfile).stem.replace(f"{basename}-", "")
    else:
        version_str = str(ncdate)
        candidates = [
            f"/u/astrid-r1b/ModelData/BedMachine/{basename}-{version_str}.nc",
            f"/home/ModelData/Antarctica/BedMachine/{basename}-{version_str}.nc",
            f"/totten_1/ModelData/Antarctica/BedMachine/{basename}-{version_str}.nc",
            f"/local/ModelData/BedMachineAntarctica/{basename}-{version_str}.nc",
            f"/Users/larour/ModelData/BedMachine/{basename}-{version_str}.nc",
            f"/Users/u7322062/Downloads/{basename}-{ncdate}.nc",
            f"./{basename}-{version_str}.nc",
        ]
        ncfile = None
        for p in candidates:
            if Path(p).exists():
                ncfile = p
                break
        if ncfile is None:
            raise FileNotFoundError(
                f"Could not find {basename}-{version_str}.nc; add its path to the list "
                "or pass a full path as ncdate."
            )

    print(f"   -- BedMachine Antarctica version: {version_str}")

    # Load axes
    with Dataset(ncfile, "r") as ds:
        xdata = np.array(ds.variables["x"][:], dtype=float)
        ydata = np.array(ds.variables["y"][:], dtype=float)

    # Window the read like MATLAB (with a small border offset)
    offset = 2
    Xr = np.asarray(X, dtype=float)
    Yr = np.asarray(Y, dtype=float)
    xmin, xmax = np.nanmin(Xr), np.nanmax(Xr)
    ymin, ymax = np.nanmin(Yr), np.nanmax(Yr)

    posx = np.where(xdata <= xmax)[0]
    if posx.size == 0:
        posx = np.array([xdata.size - 1])
    id1x = max(0, int(np.searchsorted(xdata, xmin, side="left") - 1 - offset + 1))  # mimic find >= xmin, then -offset
    id2x = min(xdata.size - 1, int(posx[-1] + offset))

    posy = np.where(ydata >= ymin)[0]
    if posy.size == 0:
        posy = np.array([ydata.size - 1])
    # in MATLAB: id1y uses first index of ydata<=ymax (so it's "left" in descending/asc mix);
    # here we assume ascending y; emulate robustly:
    id1y = max(0, int(np.searchsorted(ydata, ymin, side="left") - 1 - offset + 1))
    id2y = min(ydata.size - 1, int(np.searchsorted(ydata, ymax, side="right") - 1 + offset))

    # Safe slice ranges
    xslice = slice(id1x, id2x + 1)
    yslice = slice(id1y, id2y + 1)


    # Load the requested variable (transpose like MATLAB)
    varname = "mask" if string == "icemask" else string
    print(f"   -- BedMachine Antarctica: loading {string}")
  # --- Load the requested variable with correct axis order -----------------
    #print(f"   -- BedMachine Antarctica: loading {string}")
    with Dataset(ncfile, "r") as ds:
        if varname not in ds.variables:
            raise KeyError(f"Variable '{varname}' not found in {ncfile}")
        v = ds.variables[varname]

        # BedMachine vars can be saved as (y, x) or (x, y). Normalize to (y, x).
        dims = tuple(v.dimensions)

        # Build slices for each axis in the variable
        # We’ll map our xslice/yslice onto the correct positions.
        def _read_window_xy_as_yx(var):
            # Try to detect names 'x'/'y' in v.dimensions
            if "y" in dims and "x" in dims:
                yi = dims.index("y")
                xi = dims.index("x")
                # Construct slice tuple in var's native order
                sl = [slice(None)] * var.ndim
                sl[yi] = yslice
                sl[xi] = xslice
                arr = np.array(var[tuple(sl)], dtype=float)
                # Now ensure final array is (y, x): move axes if needed
                if yi > xi:
                    # arr axes are (..., x, y, ...). Bring y to axis 0 and x to axis 1
                    arr = np.moveaxis(arr, [yi, xi], [0, 1])
                else:
                    # arr axes are (..., y, x, ...). Bring y->0, x->1
                    arr = np.moveaxis(arr, [yi, xi], [0, 1])
                # If there were extra leading dims, squeeze them out
                arr = np.squeeze(arr)
            else:
                # Fallback: assume last 2 dims are (y, x)
                arr = np.array(var[...], dtype=float)
                arr = np.squeeze(arr)
                if arr.ndim < 2:
                    raise ValueError(f"Variable '{var.name}' is not at least 2-D.")
                arr = arr[yslice, xslice]
            return arr

        data = _read_window_xy_as_yx(v)  # data now intended to be (y, x)

    # Windowed grid vectors
    xw = xdata[xslice]
    yw = ydata[yslice]

    # MATLAB code transposed after read; we’ve already normalized to (y, x).
    # Apply 'icemask' post-processing exactly like MATLAB.
    if string == "icemask":
        data = data.copy()
        data[data == 3] = 0

    # Ensure axes are strictly ascending for interpolators. Flip data if needed.
    if xw.size > 1 and xw[1] < xw[0]:
        xw = xw[::-1]
        data = data[:, ::-1]
    if yw.size > 1 and yw[1] < yw[0]:
        yw = yw[::-1]
        data = data[::-1, :]

    # --- Interpolate (data must be (len(yw), len(xw))) ----------------------
    print(f"   -- BedMachine Antarctica: interpolating {string}")
    print(f"       -- Interpolation method: {method}")
    if string in ("mask", "source"):
        output = _interp_from_grid(xw, yw, data, Xr, Yr, method="nearest")
    else:
        output = _interp_from_grid(xw, yw, data, Xr, Yr, method=method)
    return output


def _interp_from_grid(x, y, data, X, Y, method="cubic"):
    """
    Emulates ISSM's InterpFromGrid(x,y,data,X,Y,method) with:
      - 'nearest' / 'linear' via RegularGridInterpolator (if SciPy) or FastInterp fallback
      - 'cubic' via RectBivariateSpline (SciPy). Falls back to 'linear' if SciPy unavailable.
    Inputs:
      x, y: 1D axes (ascending), data shape = (len(y), len(x))  [we've transposed to match this]
      X, Y: target arrays (same shape)
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    Z = np.asarray(data, dtype=float)
    Xp = np.asarray(X, dtype=float)
    Yp = np.asarray(Y, dtype=float)

    if method == "nearest":
        if _HAVE_SCIPY:
            rgi = RegularGridInterpolator(
                (y, x), Z, method="nearest", bounds_error=False, fill_value=np.nan
            )
            out = rgi(np.column_stack([Yp.ravel(), Xp.ravel()]))
            return out.reshape(Xp.shape)
        # Fallback: nearest via FastInterp
        return _fast_interp(x, y, Z, Xp, Yp, method="nearest")

    if method == "linear":
        if _HAVE_SCIPY:
            rgi = RegularGridInterpolator(
                (y, x), Z, method="linear", bounds_error=False, fill_value=np.nan
            )
            out = rgi(np.column_stack([Yp.ravel(), Xp.ravel()]))
            return out.reshape(Xp.shape)
        # Fallback: bilinear via FastInterp
        return _fast_interp(x, y, Z, Xp, Yp, method="bilinear")

    if method == "cubic":
        if _HAVE_SCIPY:
            # RectBivariateSpline expects strictly ascending axes
            rbs = RectBivariateSpline(y, x, Z, kx=3, ky=3)
            out = rbs.ev(Yp.ravel(), Xp.ravel()).reshape(Xp.shape)
            return out
        # No SciPy: degrade gracefully to linear
        return _fast_interp(x, y, Z, Xp, Yp, method="bilinear")

    # Unknown method -> try linear as a sane default
    if _HAVE_SCIPY:
        rgi = RegularGridInterpolator(
            (y, x), Z, method="linear", bounds_error=False, fill_value=np.nan
        )
        out = rgi(np.column_stack([Yp.ravel(), Xp.ravel()]))
        return out.reshape(Xp.shape)
    return _fast_interp(x, y, Z, Xp, Yp, method="bilinear")


def _fast_interp(x, y, Z, X, Y, method="bilinear"):
    """
    Vectorized reproduction of MATLAB FastInterp:
      - method='nearest' or 'bilinear'
    Assumes:
      - x, y are 1D ascending
      - Z has shape (len(y), len(x))  [row=y, col=x]
      - X, Y broadcastable / same shape for output
    """
    M, N = Z.shape  # rows=ny, cols=nx
    x0, y0 = x[0], y[0]
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    ndx = 1.0 / dx
    ndy = 1.0 / dy

    xi = (X - x0) * ndx
    yi = (Y - y0) * ndy

    out = np.full_like(X, np.nan, dtype=float)

    if method == "nearest":
        rxi = np.rint(xi).astype(int) + 1  # MATLAB 1-based
        ryi = np.rint(yi).astype(int) + 1
        # In-bounds mask (1..N / 1..M)
        mask = (rxi > 0) & (rxi <= N) & np.isfinite(rxi) & (ryi > 0) & (ryi <= M) & np.isfinite(ryi)
        # Convert to 0-based indices
        j = np.clip(rxi - 1, 0, N - 1)
        i = np.clip(ryi - 1, 0, M - 1)
        out[mask] = Z[i[mask], j[mask]]
        return out

    # bilinear
    fxi = np.floor(xi).astype(int) + 1  # 1-based xi
    fyi = np.floor(yi).astype(int) + 1  # 1-based yi
    dfxi = xi - (fxi - 1)
    dfyi = yi - (fyi - 1)

    mask = (
        (fxi > 0) & (fxi < N) & np.isfinite(fxi) &
        (fyi > 0) & (fyi < M) & np.isfinite(fyi)
    )
    if not np.any(mask):
        return out

    fxi0 = np.clip(fxi[mask] - 1, 0, N - 1)  # -> 0-based
    fyi0 = np.clip(fyi[mask] - 1, 0, M - 1)
    dfxi_m = dfxi[mask]
    dfyi_m = dfyi[mask]

    # Four neighbors
    i1 = fyi0
    j1 = fxi0
    i2 = fyi0
    j2 = np.clip(fxi0 + 1, 0, N - 1)
    i3 = np.clip(fyi0 + 1, 0, M - 1)
    j3 = np.clip(fxi0 + 1, 0, N - 1)
    i4 = np.clip(fyi0 + 1, 0, M - 1)
    j4 = fxi0

    z1 = Z[i1, j1]
    z2 = Z[i2, j2]
    z3 = Z[i3, j3]
    z4 = Z[i4, j4]

    out_m = (
        z1 * (1 - dfxi_m) * (1 - dfyi_m) +
        z2 * (dfxi_m) * (1 - dfyi_m) +
        z4 * (1 - dfxi_m) * (dfyi_m) +
        z3 * (dfxi_m) * (dfyi_m)
    )
    out[mask] = out_m
    return out
