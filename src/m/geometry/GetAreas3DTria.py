import numpy as np


def GetAreas3DTria(index, x, y, z, *args):
    """GETAREAS3DTRIA - compute areas of triangles with 3D coordinates

    Compute areas of triangles with 3D coordinates.

    Usage:
        areas = GetAreas3DTria(index, x, y, z)

    Examples:
        areas = GetAreas3DTria(md.mesh.elements, md.mesh.x, md.mesh.y, md.mesh.z)

    TODO:
    - Determine if *args is needed.
    """

    # Get number of elements and number of nodes
    nels = index.shape[0]
    nods = len(x)

    # Some checks
    nargs = len(args)

    # TODO: Do we really need this under Python (first 4 arguments are required)?
    # if nargs != 3 and nargs != 4:
    #     print(GetAreas3DTria.__doc__)
    #     raise Exception('GetAreas3DTria error message: bad usage')

    if len(y) != nods or (nargs == 4 and len(z) != nods):
        print(GetAreas3DTria.__doc__)
        raise Exception('GetAreas3DTria error message: x, y, and z do not have the same length')

    if np.max(index) > nods:
        print(GetAreas3DTria.__doc__)
        raise Exception('GetAreas3DTria error message: index should not have values above {}'.format(nods))

    if nargs == 4 and index.shape[1] != 3:
        print(GetAreas3DTria.__doc__)
        raise Exception('GetAreas3DTria error message: index should have 3 columns for 2d meshes')

    # Initialization
    areas = np.zeros((nels, ))
    x1 = x[index[:, 0] - 1]
    x2 = x[index[:, 1] - 1]
    x3 = x[index[:, 2] - 1]
    y1 = y[index[:, 0] - 1]
    y2 = y[index[:, 1] - 1]
    y3 = y[index[:, 2] - 1]
    z1 = z[index[:, 0] - 1]
    z2 = z[index[:, 1] - 1]
    z3 = z[index[:, 2] - 1]

    # Area of triangles with 3D coordinates
    for i in range(nels):
        m1 = np.vstack(([x1[i], x2[i], x3[i]], [y1[i], y2[i], y3[i]], [1, 1, 1]))
        m2 = np.vstack(([y1[i], y2[i], y3[i]], [z1[i], z2[i], z3[i]], [1, 1, 1]))
        m3 = np.vstack(([z1[i], z2[i], z3[i]], [x1[i], x2[i], x3[i]], [1, 1, 1]))
        areas[i] = ((np.linalg.det(m1) ** 2 + np.linalg.det(m2) ** 2 + np.linalg.det(m3) ** 2) ** 0.5) / 2 # NOTE: math.sqrt cannot be applied element-wise to a list/numpy.array

    return areas
