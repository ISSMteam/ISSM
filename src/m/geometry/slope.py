import numpy as np
from GetNodalFunctionsCoeff import GetNodalFunctionsCoeff
from project3d import project3d


def slope(md, *args):
    """
    SLOPE - compute the surface slope

    Usage:
            sx, sy, s = slope(md)
            sx, sy, s = slope(md, md.results.TransientSolution(1).Surface)
    """

    #load some variables (it is much faster if the variables are loaded from md once for all)
    if md.mesh.dimension() == 2:
        index = md.mesh.elements
        x = md.mesh.x
        y = md.mesh.y
    else:
        index = md.mesh.elements2d
        x = md.mesh.x2d
        y = md.mesh.y2d

    if len(args) == 0:
        surf = md.geometry.surface
    elif len(args) == 1:
        surf = args[0]
    else:
        raise RuntimeError("slope.py usage error")

    #%compute nodal functions coefficients N(x, y)=alpha x + beta y + gamma
    alpha, beta = GetNodalFunctionsCoeff(index, x, y)[0:2]

    summation = np.array([[1], [1], [1]])
    sx = np.dot(surf[index - 1] * alpha, summation).reshape(-1, )
    sy = np.dot(surf[index - 1] * beta, summation).reshape(-1, )

    s = np.sqrt(sx**2 + sy**2)

    if md.mesh.dimension() == 3:
        sx = project3d(md, 'vector', sx, 'type', 'element')
        sy = project3d(md, 'vector', sy, 'type', 'element')
        s = np.sqrt(sx**2 + sy**2)

    return (sx, sy, s)
