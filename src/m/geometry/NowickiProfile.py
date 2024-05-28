import numpy as np


def NowickiProfile(x):
    """
    NOWICKIPROFILE - Create profile at the transition zone based on Sophie Nowicki's thesis

    Usage:
        [b h] = NowickiProfile(x)

         - h = ice thickness
         - b = ice base
         - x = along flow coordinate
    """
    #Constant for theoretical profile
    delta = 0.1  #ratio of water density and ice density - 1
    hg = 1.  #ice thickness at grounding line
    sea = hg / (1 + delta)  #sea level
    lamda = 0.1  #ration of deviatoric stress and water pressure
    beta = 5.  #friction coefficient
    ms = 0.005  #surface accumulation rat
    mu = 5.  #viscosity
    q = 0.801  #ice mass flux

    #mesh parameters
    b = np.zeros((np.size(x), ))
    h = np.zeros((np.size(x), ))
    s = np.zeros((np.size(x), ))

    #upstream of the GL
    for i in range(int(np.size(x) / 2)):
        ss = np.roots([1, 4 * lamda * beta, 0, 0, 6 * lamda * ms * x[i]**2 + 12 * lamda * q * x[i] - hg**4 - 4 * lamda * beta * hg**3])
        for j in range(4):
            if (np.isreal(ss[j]) > 0) and (np.imag(ss[j]) == 0):
                s[i] = ss[j]
        h[i] = s[i]
        b[i] = 0.

    #downstream of the GL
    for i in range(int(np.size(x) / 2), int(np.size(x))):
        h[i] = (x[i] / (4. * (delta + 1) * q) + hg**(-2))**(-0.5)  # ice thickness for ice shelf from (3.1)
        b[i] = sea - h[i] * (1. / (1 + delta))

    return [b, h, sea]
