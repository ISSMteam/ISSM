import numpy as np


def GetNodalFunctionsCoeff(index, x, y):
    """
    GETNODELFUNCTIONSCOEFF - compute nodal functions coefficients

       Compute the coefficients alpha beta and optionaly gamma of
       2d triangular elements. For each element, the nodal function
       is defined as:
       N(x, y)=sum(i = 1:3) alpha_i * x + beta_i * y + gamma_i

       Usage:
          [alpha beta] = GetNodalFunctionsCoeff(index, x, y)
          [alpha beta gamma] = GetNodalFunctionsCoeff(index, x, y)

       Example:
          [alpha beta gamma] = GetNodalFunctionsCoeff(md.mesh.elements, md.mesh.x, md.mesh.y)
    """

    #make columns out of x and y
    x = x.reshape(-1)
    y = y.reshape(-1)

    #get nels and nods
    nels = np.size(index, axis=0)
    nods = np.size(x)

    #some checks
    if np.size(y) != nods:
        raise TypeError("GetNodalFunctionsCoeff error message: x and y do not have the same length.")
    if np.max(index) > nods:
        raise TypeError("GetNodalFunctionsCoeff error message: index should not have values above {}.".format(nods))
    if np.size(index, axis=1) != 3:
        raise TypeError("GetNodalFunctionsCoeff error message: only 2d meshes supported. index should have 3 columns.")

    #initialize output
    alpha = np.zeros((nels, 3))
    beta = np.zeros((nels, 3))

    #compute nodal functions coefficients N(x, y) = alpha x + beta y + gamma
    x1 = x[index[:, 0] - 1]
    x2 = x[index[:, 1] - 1]
    x3 = x[index[:, 2] - 1]
    y1 = y[index[:, 0] - 1]
    y2 = y[index[:, 1] - 1]
    y3 = y[index[:, 2] - 1]
    invdet = 1. / (x1 * (y2 - y3) - x2 * (y1 - y3) + x3 * (y1 - y2))

    #get alpha and beta
    alpha = np.vstack(((invdet * (y2 - y3)).reshape(-1, ), (invdet * (y3 - y1)).reshape(-1, ), (invdet * (y1 - y2)).reshape(-1, ))).T
    beta = np.vstack(((invdet * (x3 - x2)).reshape(-1, ), (invdet * (x1 - x3)).reshape(-1, ), (invdet * (x2 - x1)).reshape(-1, ))).T

    #get gamma if requested
    gamma = np.zeros((nels, 3))
    gamma = np.vstack(((invdet * (x2 * y3 - x3 * y2)).reshape(-1, ), (invdet * (y1 * x3 - y3 * x1)).reshape(-1, ), (invdet * (x1 * y2 - x2 * y1)).reshape(-1, ))).T

    return alpha, beta, gamma
