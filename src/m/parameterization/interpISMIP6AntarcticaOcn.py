
import os, re, glob
import numpy as np
from pathlib import Path
from netCDF4 import Dataset
from typing import Optional, Tuple, List

# Optional ISSM helpers
try:
    from InterpFromGridToMesh import InterpFromGridToMesh as _ISSM_InterpFromGrid
except Exception:
    _ISSM_InterpFromGrid = None

try:
    from basalforcingsismip6 import basalforcingsismip6 as _BasalForcingsISSM6
except Exception:
    _BasalForcingISSM6 = None

# Fallback interpolator
from scipy.interpolate import griddata as _griddata

def _decimal_year(year, month, day=15):
    return float(year) + (float(month) - 0.5) / 12.0
def _hostname_based_basepath() -> str:
    """Replicates your MATLAB switch on hostname; raise if unknown."""
    host = os.uname().nodename
    if host == "totten":
        return "/totten_1/ModelData/ISMIP6/Projections/AIS/Ocean_Forcing/"
    elif host == "amundsen.thayer.dartmouth.edu":
        return "/local/ModelData/ISMIP6Data/Forcings2100/Ocean/"
    elif host == "nri-085597":
        return "/Users/u7322062/Documents/ISSM/examples/mismip_reduced/DATA/ISMIP6/Ocean_Forcing/"
    raise RuntimeError(f"Machine '{host}' not supported yet; please provide your own path")


def _interp2d(
    xaxis: np.ndarray,
    yaxis: np.ndarray,
    field2d: np.ndarray,
    xq: np.ndarray,
    yq: np.ndarray,
) -> np.ndarray:
    """
    Interpolate rectilinear grid -> scattered points.

    xaxis (nx,), yaxis (ny,), field2d (ny, nx)
    """
    if _ISSM_InterpFromGrid is not None:
        # ISSM version (fast)
        # MATLAB used field' when calling, so transpose here as well.
        return _ISSM_InterpFromGrid(xaxis, yaxis, field2d.T, xq, yq, 0)

    # SciPy fallback (slower, but OK since we only do Nz interpolations)
    X, Y = np.meshgrid(xaxis, yaxis)  # (ny, nx)
    pts = np.c_[X.ravel(), Y.ravel()]
    vals = field2d.ravel()

    out = _griddata(pts, vals, np.c_[xq, yq], method="linear")
    if np.isnan(out).any():
        out_near = _griddata(pts, vals, np.c_[xq, yq], method="nearest")
        out = np.where(np.isnan(out), out_near, out)
    return out


def _months_1995_2014() -> List[Tuple[int, int]]:
    months: List[Tuple[int, int]] = []
    for y in range(1995, 2015):  # inclusive 1995–2014
        for m in range(1, 13):
            months.append((y, m))
    return months  # length 240

from model import *
# --------- main function ---------
def interpISMIP6AntarcticaOcn(
    md: Optional[object] = None,
    model_name: str = "",
    root: Optional[str] = None,
    filename: str = "obs_thermal_forcing_1995-2017_8km_x_60m.nc",
):
    """
    Load observational TF from a single 3-D file thermal_forcing(z,y,x),
    build a monthly time axis 1995–2014, and replicate the 3-D field
    for each month (constant in time).

    Returns an object suitable for md.basalforcings (basalforcingsismip6).
    """
    class _BF(object):
        """Dummy basalforcings object for when basalforcingsismip6 Python wrapper
        is not available. Must expose a 'checkconsistency' attribute so that
        loadvars/ismodelselfconsistent don't crash.
        """
        def __init__(self):
            # Simple no-op consistency check
            def _checkconsistency(md, solution, analyses):
                return
            self.checkconsistency = _checkconsistency
    # --- Case 1: called by loadmodel with no args ---------------------------
     # --- Called by loadmodel/loadvars with no args ---------------------------
    if md is None:
        if _BasalForcingsISSM6 is not None:
            return _BasalForcingsISSM6()
        else:
            return _BF()   # <-- important




    # ----- locate file -----
    base = _hostname_based_basepath() if root is None else root
    candidate_paths = [
        Path(base) / model_name / "1995-2017" / filename,
        Path(base) / model_name / "1995-2100" / filename,
        Path(base) / model_name / filename,
        Path(filename),
    ]
    ncpath = next((p for p in candidate_paths if p.exists()), None)
    if ncpath is None:
        raise FileNotFoundError(f"Could not find {filename} under {base}/{model_name}")

    ncpath = str(ncpath)

    # ----- read axes + 3D TF (no time) -----
    with Dataset(ncpath, "r") as nc:
        x_n = np.array(nc.variables["x"][:], dtype=float)
        y_n = np.array(nc.variables["y"][:], dtype=float)
        z_data = np.array(nc.variables["z"][:], dtype=float)
        tf3d = np.array(nc.variables["thermal_forcing"][:], dtype=float)  # (z,y,x)

    Nz, Ny, Nx = tf3d.shape

    # ----- build time: monthly 1995–2014 -----
    months = _months_1995_2014()
    Nt = len(months)  # 240

        # ----- build time: ONLY ONE TIME since field is constant -----
    # pick an arbitrary reference time, e.g. mid-1995
    time = np.array([_decimal_year(1995, 7)], dtype=float)  # shape (1,)

    Nnode = md.mesh.numberofvertices
    tf: List[np.ndarray] = [None] * Nz

    print(
        f"   == Interpolating observational TF (constant in time): "
        f"Nz={Nz}, Nt={time.size}, nodes={Nnode}"
    )

    # ----- interpolate per depth -----
    for iz in range(Nz):
        print(f"   == Interpolating over depth {iz+1}/{Nz}")
        Fyx = tf3d[iz, :, :]  # (y,x)

        # Interpolate ONCE at this depth
        vals = _interp2d(x_n, y_n, Fyx, md.mesh.x, md.mesh.y).astype(np.float32)  # (Nnode,)

        # Single time slice: (Nnode, 1)
        temp_matrix = vals.reshape(-1, 1)

        # Append a single time row at bottom: (Nnode+1, 1)
        tf_depth = np.vstack([temp_matrix, time.reshape(1, -1)])  # (Nnode+1, 1)
        tf[iz] = tf_depth

    # time = np.array([_decimal_year(y, m) for (y, m) in months], dtype=float)

    # Nnode = md.mesh.numberofvertices
    # tf: List[np.ndarray] = [None] * Nz

    # print(
    #     f"   == Interpolating observational TF (constant in time) 1995–2014: "
    #     f"Nz={Nz}, Nt={Nt}, nodes={Nnode}"
    # )

    # # ----- interpolate per depth -----
    # for iz in range(Nz):
    #     print(f"   == Interpolating over depth {iz+1}/{Nz}")
    #     Fyx = tf3d[iz, :, :]  # (y,x)

    #     # Interpolate ONCE at this depth (field is time-invariant)
    #     vals = _interp2d(x_n, y_n, Fyx, md.mesh.x, md.mesh.y)  # (Nnode,)

    #     # Replicate this column for all months: (Nnode, Nt)
    #     temp_matrix = np.tile(vals.reshape(-1, 1), (1, Nt))

    #     # Append time row at bottom, like MATLAB [values ; time]
    #     tf_depth = np.vstack([temp_matrix, time.reshape(1, -1)])
    #     tf[iz] = tf_depth

    # ---------- Basin / ΔT / γ0 handling (same as before) ----------
    path = _hostname_based_basepath() if root is None else root
    deltatnc_median = str(
        Path(path)
        / "parameterizations"
        / "coeff_gamma0_DeltaT_quadratic_non_local_median.nc"
    )
    basin_datanc = str(Path(path) / "imbie2" / "imbie2_basin_numbers_8km.nc")

    with Dataset(deltatnc_median, "r") as nc:
        deltaT_median = np.array(nc["deltaT_basin"][:], dtype=float)
        gamma0_median = np.array(nc["gamma0"][:], dtype=float)
    with Dataset(basin_datanc, "r") as nc:
        basinid_grid = np.array(nc["basinNumber"][:], dtype=float)

    unique_ids = np.unique(basinid_grid[~np.isnan(basinid_grid)]).astype(int)
    num_basins = unique_ids.size
    deltat_median = np.empty(num_basins, dtype=float)
    for i, bid in enumerate(unique_ids):
        vals = deltaT_median[basinid_grid == bid]
        deltat_median[i] = vals.flat[0] if vals.size else np.nan

    # element-center basin id (handle 0- vs 1-based elements)
    if getattr(md.mesh, "elements", None) is None:
        raise RuntimeError("md.mesh.elements is required to compute basin IDs")

    one_based = np.min(md.mesh.elements) == 1
    elems = md.mesh.elements - 1 if one_based else md.mesh.elements
    x_el = np.mean(md.mesh.x[elems], axis=1)
    y_el = np.mean(md.mesh.y[elems], axis=1)

    # ensure basin grid shape matches (ny,nx)
    if basinid_grid.shape == (Ny, Nx):
        basin_yx = basinid_grid
    elif basinid_grid.shape == (Nx, Ny):
        basin_yx = basinid_grid.T
    else:
        raise RuntimeError("Unexpected basin grid shape vs (ny,nx) from obs TF file")

    Xg, Yg = np.meshgrid(x_n, y_n)
    b_near = _griddata(
        np.c_[Xg.ravel(), Yg.ravel()],
        basin_yx.ravel(),
        np.c_[x_el, y_el],
        method="nearest",
    )
    basinid = b_near.astype(int) + 1  # +1 to match MATLAB behavior

    # # ----- pack basalforcings object -----
    # if _BasalForcingsISSM6 is not None:
    #     # If md already has basalforcings, copy generic fields, otherwise new
    #     if getattr(md, "basalforcings", None) is not None:
    #         bf = _BasalForcingsISSM6(md.basalforcings)
    #     else:
    #         bf = _BasalForcingsISSM6()
    # else:
    #     bf = _BF()   # <-- use dummy with checkconsistency
    bf = md.basalforcings

    bf.basin_id = basinid
    bf.num_basins = int(num_basins)
    bf.delta_t = deltat_median
    bf.tf_depths = z_data.astype(float)
    bf.gamma_0 = gamma0_median.astype(float)
    bf.tf = tf


    print(
        f"Info: observational TF assumed constant in time; "
        f"monthly forcings constructed for 1995–2014 (Nt={Nt}) from {Path(ncpath).name}"
    )
    return md
