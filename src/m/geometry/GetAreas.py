import numpy as np


def GetAreas(index, x, y, z=np.array([])):
    """GetAreas - compute areas or volumes of elements

    compute areas of triangular elements or volumes of pentahedrons

    Usage:
        areas = GetAreas(index, x, y)
        volumes = GetAreas(index, x, y, z)

    Examples:
        areas = GetAreas(md.mesh.elements, md.mesh.x, md.mesh.y)
        volumes = GetAreas(md.mesh.elements, md.mesh.x, md.mesh.y, md.z)
    """

    # Get number of elements and number of nodes
    nels = np.shape(index)[0]
    nods = np.shape(x)[0]

    # Some checks
    if (np.shape(y)[0] != nods) or (z.size > 0 and np.shape(z)[0] != nods):
        raise TypeError('GetAreas error message: x, y and z do not have the same length.')
    if np.max(index) > nods:
        raise TypeError('GetAreas error message: index should not have values above {}.'.format(nods))
    if z.size == 0 and np.shape(index)[1] != 3:
        raise TypeError('GetAreas error message: index should have 3 columns for 2d meshes.')
    if z.size > 0 and np.shape(index)[1] != 6:
        raise TypeError('GetAreas error message: index should have 6 columns for 3d meshes.')

    # Initialization
    areas = np.zeros(nels)
    x1 = x[index[:, 0] - 1]
    x2 = x[index[:, 1] - 1]
    x3 = x[index[:, 2] - 1]
    y1 = y[index[:, 0] - 1]
    y2 = y[index[:, 1] - 1]
    y3 = y[index[:, 2] - 1]

    # Compute the volume of each element
    if z.size == 0:
        # Compute the surface of the triangle
        areas = (0.5 * ((x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1)))
    else:
        # V = area(triangle) * 1/3(z1 + z2 + z3)
        thickness = np.mean(z[index[:, 3:6] - 1]) - np.mean(z[index[:, 0:3] - 1])
        areas = (0.5 * ((x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1))) * thickness

    return areas
