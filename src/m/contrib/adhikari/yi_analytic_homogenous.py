import numpy as np

from math import pi

def love_analytic(param):
    # param = required parameters
    # param.rho = earth density [kg/m3]
    # param.mu = shear modulus [Pa]
    # param.G = universal gravitational constant [m3 kg-1 s-2] 
    # param.radius = radius of the surface and internal layers of the planet. [m]
    # param.g0 = acceleration due to gravity on the planetary surface. [m s-2]
    # param.source = model forcing type [9 volumentaric potential, 11 surface loads]
    # param.degree = SH degree

    rho = param['rho']
    mu = param['mu']
    G = param['G']
    a = np.max(param['radius'])
    g0 = param['g0']
    r = param['radius']
    source = param['source']
    n = param['degree']

    C1 = 0
    C2 = 0
    C6 = 0
    if (source == 9):
        C3 = -1 / 2 * n * (n + 1) * rho / ((2 * n ** 2 + 4 * n + 3) * mu + n * rho * g0 * a) * a ** (-n)
        C4 = 1 / 2 / (n - 1) * n ** 2 * (n + 2) * rho / ((2 * n ** 2 + 4 * n + 3) * mu + n * rho * g0 * a) * a ** (-n + 2)
        C5 = (1 + 3 / 2 / (n - 1) * n * rho * g0 * a / ((2 * n ** 2 + 4 * n + 3) * mu + n * rho * g0 * a)) * a ** (-n)
    elif (source == 11):
        denom = a ** n * (4 * G * pi * a ** 2 * n * rho ** 2 + 3 * mu * (2 * n ** 2 + 4 * n + 3))
        C3 = (n * rho * (n ** 2 - 1)) / denom
        C4 = -(a ** 2 * n ** 2 * rho * (n + 2)) / denom
        C5 = (3 * mu * (2 * n ** 2 + 4 * n + 3)) / denom

    y1 = np.divide(C1, np.power(r, n)) + np.divide(C2, np.power(r, n + 2)) + C3 * np.power(r, n + 1) + C4 * np.power(r, n - 1)
    y2 = 2 * mu * ((-n ** 2 - 3 * n + 1) / (n + 1) * C1 * np.power(r, -n - 1) - (n + 2) * C2 * np.power(r, -n - 3) + (n ** 2 - n - 3) / n * C3 * np.power(r, n) + (n - 1) * C4 * np.power(r, n - 2)) + 4 / 3 * pi * G * rho ** 2 * (C1 * np.power(r, -n + 1) + C2 * np.power(r, -n - 1) + C3 * np.power(r, n + 2) + C4 * np.power(r, n)) - rho * C5 * np.power(r, n) - rho * C6 * np.power(r, -n - 1)
    y3 = -(n - 2) / (n * (n + 1)) * C1 * np.power(r, -n) - 1 / (n + 1) * C2 * np.power(r, -n - 2) + (n + 3) / (n * (n + 1)) * C3 * np.power(r, n + 1) + C4 / n * np.power(r, n - 1)
    y4 = 2 * mu * ((n - 1) / n * C1 * np.power(r, -n - 1) + (n + 2) / (n + 1) * C2 * np.power(r, -n - 3) + (n + 2) / (n + 1) * C3 * np.power(r, n)       + (n - 1) / n * C4 * np.power(r, n - 2))
    y5 = C5 * np.power(r, n) + C6 * np.power(r, -n - 1)
    y6 = n * C5 * np.power(r, n - 1) - (n + 1) * C6 * np.power(r, -n - 2) - 4 * pi * G * rho * (C1 * np.power(r, -n) + C2 * np.power(r, -n - 2) + C3 * np.power(r, n + 1) + C4 * np.power(r, n - 1))

    return np.squeeze(np.stack((y1, y2, y3, y4, y5, y6), axis=1))
