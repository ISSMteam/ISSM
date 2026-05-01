import numpy as np
from GetNodalFunctionsCoeff import GetNodalFunctionsCoeff
from project3d import project3d

def slope(md, *args):
    """
    SLOPE - compute the gradient of any field

    Usage:
            dfdx, dfdy, ds = slope(md)
            dfdx, dfdy, ds = slope(md, md.results.TransientSolution(1).Surface)
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
        f = md.geometry.surface
    elif len(args) == 1:
        f = args[0]
    else:
        raise RuntimeError("slope.py usage error")

    #%compute nodal functions coefficients N(x, y)=alpha x + beta y + gamma
    alpha, beta = GetNodalFunctionsCoeff(index, x, y)[0:2]

    summation = np.array([[1], [1], [1]])
    dfdx = np.dot(f[index - 1] * alpha, summation).reshape(-1, )
    dfdy = np.dot(f[index - 1] * beta, summation).reshape(-1, )
    if md.mesh.dimension() == 3:
        dfdx = project3d(md, 'vector', dfdx, 'type', 'element')
        dfdy = project3d(md, 'vector', dfdy, 'type', 'element')

    #Compute magnitude
    ds = np.sqrt(dfdx**2 + dfdy**2)

    return (dfdx, dfdy, ds)
