import numpy as np


def GetAreas(index, x, y, z=np.array([])):
    """
    GETAREAS - compute areas or volumes of elements

       compute areas of triangular elements or volumes
       of pentahedrons

       Usage:
          areas  =GetAreas(index, x, y)
          volumes = GetAreas(index, x, y, z)

       Examples:
          areas  =GetAreas(md.mesh.elements, md.mesh.x, md.mesh.y)
          volumes = GetAreas(md.mesh.elements, md.mesh.x, md.mesh.y, md.z)
    """

    #get number of elements and number of nodes
    nels = np.size(index, axis=0)
    nods = np.size(x)

    #some checks
    if np.size(y) != nods or (z and np.size(z) != nods):
        raise TypeError("GetAreas error message: x, y and z do not have the same length.")
    if np.max(index) > nods:
        raise TypeError("GetAreas error message: index should not have values above %d." % nods)
    if (not z and np.size(index, axis=1) != 3):
        raise TypeError("GetAreas error message: index should have 3 columns for 2d meshes.")
    if (z and np.size(index, axis=1) != 6):
        raise TypeError("GetAreas error message: index should have 6 columns for 3d meshes.")

    #initialization
    areas = np.zeros(nels)
    x1 = x[index[:, 0] - 1]
    x2 = x[index[:, 1] - 1]
    x3 = x[index[:, 2] - 1]
    y1 = y[index[:, 0] - 1]
    y2 = y[index[:, 1] - 1]
    y3 = y[index[:, 2] - 1]

    #compute the volume of each element
    if not z:
        #compute the surface of the triangle
        areas = (0.5 * ((x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1)))
    else:
        #V = area(triangle) * 1 / 3(z1 + z2 + z3)
        thickness = np.mean(z[index[:, 3:6] - 1]) - np.mean(z[index[:, 0:3] - 1])
        areas = (0.5 * ((x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1))) * thickness

    return areas
