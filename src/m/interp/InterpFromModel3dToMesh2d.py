import numpy as np

# Expect these ISSM helpers to be importable in your environment
from DepthAverage import DepthAverage
from InterpFromMeshToMesh2d import InterpFromMeshToMesh2d
from InterpFromMeshToMesh3d import InterpFromMeshToMesh3d


def InterpFromModel3dToMesh2d(md, data, x_prime, y_prime, sigma, default_value):
    """
    Interpolate from a 3D hexahedron mesh (model `md`) onto a list of 2D points.

    Parameters
    ----------
    md : ISSM model (3D)
        Model holding mesh/geometry.
    data : array-like
        Field defined on the 3D mesh (per element or per vertex).
        Length must equal md.mesh.numberofelements or md.mesh.numberofvertices.
    x_prime, y_prime : array-like
        Coordinates of the 2D target points (same shape).
    sigma : float or np.nan
        Scaled vertical coordinate in [0, 1] (0 = base, 1 = surface).
        If NaN, do vertical average first, then 2D interpolation.
    default_value : float
        Fill value when interpolation falls in holes.

    Returns
    -------
    data_prime : np.ndarray
        Interpolated values at (x_prime, y_prime).
    """
    # ---- basic checks --------------------------------------------------------
    if not hasattr(md.mesh, "z"):
        raise ValueError("Model should be 3D (mesh is missing 'z').")

    x_prime = np.asarray(x_prime, dtype=float)
    y_prime = np.asarray(y_prime, dtype=float)
    if x_prime.shape != y_prime.shape:
        raise ValueError("x_prime and y_prime must have the same shape.")

    data = np.asarray(data)
    if (data.size != getattr(md.mesh, "numberofelements", None) and
        data.size != getattr(md.mesh, "numberofvertices", None)):
        raise ValueError("Data length must match number of elements or number of vertices.")

    # sigma validation: allow NaN (vertical average) or 0..1
    if not (np.isnan(sigma) or (0.0 <= float(sigma) <= 1.0)):
        raise ValueError("sigma must be between 0 and 1, or NaN for vertical average.")

    # ---- vertical average path ----------------------------------------------
    if np.isnan(sigma):
        # Average vertically onto the 2D mesh, then 2D interpolate to points
        averaged_data = DepthAverage(md, data)
        return InterpFromMeshToMesh2d(
            md.mesh.elements2d,
            md.mesh.x2d,
            md.mesh.y2d,
            averaged_data,
            x_prime,
            y_prime,
            default_value
        )

    # ---- 3D sampling path (sigma in [0,1]) ----------------------------------
    # Scaled vertical coordinate at 3D vertices
    thickness = np.asarray(md.geometry.thickness, dtype=float)
    base = np.asarray(md.geometry.base, dtype=float)
    z = np.asarray(md.mesh.z, dtype=float)

    # Avoid division by zero where thickness is 0 (rare/invalid); follow MATLAB’s strictness
    if np.any(thickness == 0):
        raise ValueError("Encountered zero thickness; check model geometry before interpolation.")

    alpha = (z - base) / thickness  # expected in [0, 1]

    # MATLAB code errors if alpha goes out of bounds; do the same check here
    if np.any((alpha < 0.0) | (alpha > 1.0)):
        raise ValueError("Computed alpha outside [0,1]. Check model geometry (base/thickness).")

    # Build target z along sigma plane; n points = len(x_prime)
    z_prime = np.full(x_prime.shape, float(sigma), dtype=float)

    # Keep the 2D plane strictly inside the 3D volume like MATLAB (eps shift)
    eps = np.finfo(float).eps
    if float(sigma) == 0.0:
        z_prime = z_prime + eps
    elif float(sigma) == 1.0:
        z_prime = z_prime - eps

    # Call 3D interpolation (elements, x, y, alpha, data, xq, yq, zq, default)
    return InterpFromMeshToMesh3d(
        md.mesh.elements,
        md.mesh.x,
        md.mesh.y,
        alpha,
        data,
        x_prime,
        y_prime,
        z_prime,
        default_value
    )