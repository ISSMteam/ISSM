import numpy as np
from slope import slope

def drivingstress(md):
    """Driving stress = rho_ice * g * H * slope, evaluated on elements.
    Returns px, py, pmag (each of length numberofelements)."""
    sx, sy, _ = slope(md)

    index = np.asarray(md.mesh.elements).astype(int) - 1
    thickness = np.asarray(md.geometry.thickness).flatten()

    # average thickness over elements
    thickness_bar = np.mean(thickness[index],axis=1)

    px = -md.materials.rho_ice * md.constants.g * thickness_bar * sx
    py = -md.materials.rho_ice * md.constants.g * thickness_bar * sy
    pmag = np.sqrt(px**2 + py**2)
    return px, py, pmag