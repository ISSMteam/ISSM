import numpy as np


def cuffey(temperature):
    """CUFFEY - calculates ice rigidity as a function of temperature

    rigidity (in s^(1/3)Pa) is the flow law parameter in the flow law 
    sigma = B * e(1/3) (Cuffey and Paterson, p75).

    temperature is in Kelvin degrees

    Usage:
        rigidity = cuffey(temperature)
    """

    if np.any(temperature < 0.0):
        raise RuntimeError('input temperature should be in Kelvin (positive)')

    if np.ndim(temperature) == 2:
        #T = temperature.reshape(-1, ) - 273.15
        T = temperature.flatten() - 273.15
    elif isinstance(temperature, float) or isinstance(temperature, int):
        T = np.array([temperature]) - 273.15
    else:
        T = temperature - 273.15

    rigidity = np.zeros_like(T)
    pos = np.nonzero(T <= -45)
    if len(pos):
        rigidity[pos] = pow(10, 8) * (-0.000396645116301 * (T[pos] + 50)**3 + 0.013345579471334 * (T[pos] + 50)**2 - 0.356868703259105 * (T[pos] + 50) + 7.272363035371383)
    pos = np.nonzero(np.logical_and(-45 <= T, T < -40))
    if len(pos):
        rigidity[pos] = pow(10, 8) * (-0.000396645116301 * (T[pos] + 45)**3 + 0.007395902726819 * (T[pos] + 45)**2 - 0.253161292268336 * (T[pos] + 45) + 5.772078366321591)
    pos = np.nonzero(np.logical_and(-40 <= T, T < -35))
    if len(pos):
        rigidity[pos] = pow(10, 8) * (0.000408322072669 * (T[pos] + 40)**3 + 0.001446225982305 * (T[pos] + 40)**2 - 0.208950648722716 * (T[pos] + 40) + 4.641588833612773)
    pos = np.nonzero(np.logical_and(-35 <= T, T < -30))
    if len(pos):
        rigidity[pos] = pow(10, 8) * (-0.000423888728124 * (T[pos] + 35)**3 + 0.007571057072334 * (T[pos] + 35)**2 - 0.163864233449525 * (T[pos] + 35) + 3.684031498640382)
    pos = np.nonzero(np.logical_and(-30 <= T, T < -25))
    if len(pos):
        rigidity[pos] = pow(10, 8) * (0.000147154327025 * (T[pos] + 30)**3 + 0.001212726150476 * (T[pos] + 30)**2 - 0.119945317335478 * (T[pos] + 30) + 3.001000667185614)
    pos = np.nonzero(np.logical_and(-25 <= T, T < -20))
    if len(pos):
        rigidity[pos] = pow(10, 8) * (-0.000193435838672 * (T[pos] + 25)**3 + 0.003420041055847 * (T[pos] + 25)**2 - 0.096781481303861 * (T[pos] + 25) + 2.449986525148220)
    pos = np.nonzero(np.logical_and(-20 <= T, T < -15))
    if len(pos):
        rigidity[pos] = pow(10, 8) * (0.000219771255067 * (T[pos] + 20)**3 + 0.000518503475772 * (T[pos] + 20)**2 - 0.077088758645767 * (T[pos] + 20) + 2.027400665191131)
    pos = np.nonzero(np.logical_and(-15 <= T, T < -10))
    if len(pos):
        rigidity[pos] = pow(10, 8) * (-0.000653438900191 * (T[pos] + 15)**3 + 0.003815072301777 * (T[pos] + 15)**2 - 0.055420879758021 * (T[pos] + 15) + 1.682390865739973)
    pos = np.nonzero(np.logical_and(-10 <= T, T < -5))
    if len(pos):
        rigidity[pos] = pow(10, 8) * (0.000692439419762 * (T[pos] + 10)**3 - 0.005986511201093 * (T[pos] + 10)**2 - 0.066278074254598 * (T[pos] + 10) + 1.418983411970382)
    pos = np.nonzero(np.logical_and(-5 <= T, T < -2))
    if len(pos):
        rigidity[pos] = pow(10, 8) * (-0.000132282004110 * (T[pos] + 5)**3 + 0.004400080095332 * (T[pos] + 5)**2 - 0.074210229783403 * (T[pos] + 5) + 1.024485188140279)
    pos = np.nonzero(-2 <= T)
    if len(pos):
        rigidity[pos] = pow(10, 8) * (-0.000132282004110 * (T[pos] + 2)**3 + 0.003209542058346 * (T[pos] + 2)**2 - 0.051381363322371 * (T[pos] + 2) + 0.837883605537096)

    # Now make sure that rigidity is positive
    pos = np.nonzero(rigidity < 0)
    rigidity[pos] = pow(1, 6)

    return rigidity
