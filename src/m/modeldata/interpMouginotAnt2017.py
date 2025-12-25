import os
import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import RegularGridInterpolator

def interp_mouginot_ant2017(X, Y, *, speed=False, method="nearest", fill_value=np.nan):
    """
    Interpolate Mouginot et al. (Antarctica) velocities onto points (X, Y).

    Parameters
    ----------
    X, Y : array_like
        Target coordinates (same shape). These should be in the dataset's x/y
        coordinate system (EPSG:3031 for the common NSIDC grid).
    speed : bool, default False
        If True, return speed (sqrt(VX^2 + VY^2)) as a single array.
        If False, return (VX, VY).
    method : {"linear","nearest"}
        Interpolation method.
    fill_value : float
        Value for points outside the data domain.

    Returns
    -------
    VX, VY : ndarray (if speed=False)
        Interpolated components matching the shape of X.
    or
    S : ndarray (if speed=True)
        Interpolated speed.
    """
    # --- Choose file by hostname (like MATLAB switch) -------------------------
    host = os.uname().nodename
    if host == "ronne":
        nc_path = "/home/ModelData/Antarctica/MouginotVel/vel_nsidc.CF16_2.nc"
    elif host == "totten":
        nc_path = "/totten_1/ModelData/Antarctica/MouginotVel/vel_nsidc.CF16_2.nc"
    elif host == "amundsen.thayer.dartmouth.edu":
        nc_path = "/local/ModelData/AntarcticVelocity/v_mix.v13Mar2019.nc"
    elif host == "nri-085597":
        nc_path = "/Users/u7322062/Downloads/antarctica_ice_velocity_450m_v2.nc"
    else:
        raise RuntimeError(f"hostname not supported yet: {host!r}")

    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    if X.shape != Y.shape:
        raise ValueError("X and Y must have the same shape.")

    offset = 2  # match MATLAB

    with Dataset(nc_path, "r") as ds:
        xdata = np.array(ds.variables["x"][:], dtype=float)  # (nx,)
        ydata = np.array(ds.variables["y"][:], dtype=float)  # (ny,)

        # Determine if axes are ascending; if not, remember to flip slices.
        x_asc = np.all(np.diff(xdata) > 0)
        y_asc = np.all(np.diff(ydata) > 0)

        xmin, xmax = np.nanmin(X), np.nanmax(X)
        ymin, ymax = np.nanmin(Y), np.nanmax(Y)

        # Helper to get [i1, i2] inclusive bounds with padding for a 1D axis
        def _bounds(axis, vmin, vmax, ascending):
            n = axis.size
            if ascending:
                i1 = max(0, int(np.searchsorted(axis, vmin, side="left")) - offset)
                i2 = min(n - 1, int(np.searchsorted(axis, vmax, side="right") - 1 + offset))
            else:
                # For descending arrays, invert the logic
                # Use -axis to leverage searchsorted on ascending data
                neg = -axis
                i1 = max(0, int(np.searchsorted(neg, -vmax, side="left")) - offset)
                i2 = min(n - 1, int(np.searchsorted(neg, -vmin, side="right") - 1 + offset))
            if i1 > i2:
                raise ValueError("Requested window is outside data domain.")
            return i1, i2

        id1x, id2x = _bounds(xdata, xmin, xmax, x_asc)
        id1y, id2y = _bounds(ydata, ymin, ymax, y_asc)

        # Slice coordinate vectors
        xs = xdata[id1x:id2x + 1]
        ys = ydata[id1y:id2y + 1]

        # Read VX/VY slices; file variables are (y, x) or (x, y) depending on product.
        # The MATLAB code transposes after read, implying file layout was (x,y) -> they wanted (y,x).
        # We'll read in [y, x] order and ensure arrays end up as (ny, nx) with ys, xs.
        # netCDF4 slicing is [y, x] if variable dims are ('y','x'); we'll try that first.
        def _read_var(name):
            var = ds.variables[name]
            # Guess dim order
            dims = var.dimensions  # tuple of names
            if len(dims) != 2:
                raise ValueError(f"Variable {name} has unexpected dims: {dims}")
            # Build slices per order
            if dims[0].lower().startswith("y") and dims[1].lower().startswith("x"):
                arr = np.array(var[id1y:id2y + 1, id1x:id2x + 1], dtype=float)
            elif dims[0].lower().startswith("x") and dims[1].lower().startswith("y"):
                arr = np.array(var[id1x:id2x + 1, id1y:id2y + 1], dtype=float).T
            else:
                # Fallback: assume (y,x)
                arr = np.array(var[id1y:id2y + 1, id1x:id2x + 1], dtype=float)
            return arr

        vx_slice = _read_var("VX")
        vy_slice = _read_var("VY")

        # Make sure axes are ascending for RegularGridInterpolator
        if not x_asc:
            xs = xs[::-1]
            vx_slice = vx_slice[:, ::-1]
            vy_slice = vy_slice[:, ::-1]
        if not y_asc:
            ys = ys[::-1]
            vx_slice = vx_slice[::-1, :]
            vy_slice = vy_slice[::-1, :]

    # Build interpolators on (ys, xs) -> data[ny, nx]
    rgi_vx = RegularGridInterpolator((ys, xs), vx_slice, method=method,
                                     bounds_error=False, fill_value=fill_value)
    rgi_vy = RegularGridInterpolator((ys, xs), vy_slice, method=method,
                                     bounds_error=False, fill_value=fill_value)

    pts = np.column_stack([Y.ravel(), X.ravel()])  # note order: (y, x)
    vx = rgi_vx(pts).reshape(X.shape)
    vy = rgi_vy(pts).reshape(X.shape)

    if speed:
        return np.hypot(vx, vy)
    return vx, vy
