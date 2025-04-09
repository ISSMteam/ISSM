import numpy as np


def paterson(temperature):
    """PATERSON - figure out the rigidity of ice for a given temperature

    rigidity (in s^(1/3)Pa) is the flow law paramter in the flow law 
    sigma = B * e(1/3) (Paterson, p97).

    temperature is in Kelvin degrees

    Usage:
        rigidity = paterson(temperature)
    """

    print('paterson is outdated, please consider using cuffey instead')

    if np.any(temperature < 0.0):
        raise RuntimeError('input temperature should be in Kelvin (positive)')

    if np.ndim(temperature) == 2:
        #T = temperature.reshape(-1, ) - 273.15
        T = temperature.flatten() - 273.15
    elif isinstance(temperature, float) or isinstance(temperature, int):
        T = np.array([temperature]) - 273.15
    else:
        T = temperature - 273.15

    # The routine below is equivalent to:

    # n = 3; T = temperature-273
    # %From paterson,
    # Temp = [0; -2; -5; -10; -15; -20; -25; -30; -35; -40; -45; -50]
    # A = [6.8 * 1.0e-15;2.4 * 1.0e-15;1.6 * 1.0e-15;4.9 * 1.0e-16;2.9 * 1.0e-16;1.7 * 1.0e-16;9.4 *
    # 1.0e-17;5.1 * 1.0e-17;2.7 * 1.0e-17;1.4 * 1.0e-17;7.3 * 1.0e-18;3.6 * 1.0e-18];;%s - 1(kPa - 3)
    # %Convert into rigidity B
    # B = A.^(-1 / n) * 1.0e3; %s^(1 / 3)Pa
    # %Now, do a cubic fit between Temp and B:
    # fittedmodel = fit(Temp, B, 'cubicspline')
    # rigidity = fittedmodel(temperature)

    rigidity = np.zeros_like(T)
    pos1 = np.nonzero(T <= -45)[0]
    if len(pos1):
        rigidity[pos1] = pow(10, 8) * (-0.000292866376675 * (T[pos1] + 50)**3 + 0.011672640664130 * (T[pos1] + 50)**2 - 0.325004442485481 * (T[pos1] + 50) + 6.524779401948101)
    pos2 = np.nonzero(np.logical_and(-45 <= T, T < -40))[0]
    if len(pos2):
        rigidity[pos2] = pow(10, 8) * (-0.000292866376675 * (T[pos2] + 45)**3 + 0.007279645014004 * (T[pos2] + 45)**2 - 0.230243014094813 * (T[pos2] + 45) + 5.154964909039554)
    pos3 = np.nonzero(np.logical_and(-40 <= T, T < -35))[0]
    if len(pos3):
        rigidity[pos3] = pow(10, 8) * (0.000072737147457 * (T[pos3] + 40)**3 + 0.002886649363879 * (T[pos3] + 40)**2 - 0.179411542205399 * (T[pos3] + 40) + 4.149132666831214)
    pos4 = np.nonzero(np.logical_and(-35 <= T, T < -30))[0]
    if len(pos4):
        rigidity[pos4] = pow(10, 8) * (-0.000086144770023 * (T[pos4] + 35)**3 + 0.003977706575736 * (T[pos4] + 35)**2 - 0.145089762507325 * (T[pos4] + 35) + 3.333333333333331)
    pos5 = np.nonzero(np.logical_and(-30 <= T, T < -25))[0]
    if len(pos5):
        rigidity[pos5] = pow(10, 8) * (-0.000043984685769 * (T[pos5] + 30)**3 + 0.002685535025386 * (T[pos5] + 30)**2 - 0.111773554501713 * (T[pos5] + 30) + 2.696559088937191)
    pos6 = np.nonzero(np.logical_and(-25 <= T, T < -20))[0]
    if len(pos6):
        rigidity[pos6] = pow(10, 8) * (-0.000029799523463 * (T[pos6] + 25)**3 + 0.002025764738854 * (T[pos6] + 25)**2 - 0.088217055680511 * (T[pos6] + 25) + 2.199331606342181)
    pos7 = np.nonzero(np.logical_and(-20 <= T, T < -15))[0]
    if len(pos7):
        rigidity[pos7] = pow(10, 8) * (0.000136920904777 * (T[pos7] + 20)**3 + 0.001578771886910 * (T[pos7] + 20)**2 - 0.070194372551690 * (T[pos7] + 20) + 1.805165505978111)
    pos8 = np.nonzero(np.logical_and(-15 <= T, T < -10))[0]
    if len(pos8):
        rigidity[pos8] = pow(10, 8) * (-0.000899763781026 * (T[pos8] + 15)**3 + 0.003632585458564 * (T[pos8] + 15)**2 - 0.044137585824322 * (T[pos8] + 15) + 1.510778053489523)
    pos9 = np.nonzero(np.logical_and(-10 <= T, T < -5))[0]
    if len(pos9):
        rigidity[pos9] = pow(10, 8) * (0.001676964325070 * (T[pos9] + 10)**3 - 0.009863871256831 * (T[pos9] + 10)**2 - 0.075294014815659 * (T[pos9] + 10) + 1.268434288203714)
    pos10 = np.nonzero(np.logical_and(-5 <= T, T < -2))[0]
    if len(pos10):
        rigidity[pos10] = pow(10, 8) * (-0.003748937622487 * (T[pos10] + 5)**3 + 0.015290593619213 * (T[pos10] + 5)**2 - 0.048160403003748 * (T[pos10] + 5) + 0.854987973338348)
    pos11 = np.nonzero(-2 <= T)[0]
    if len(pos11):
        rigidity[pos11] = pow(10, 8) * (-0.003748937622488 * (T[pos11] + 2)**3 - 0.018449844983174 * (T[pos11] + 2)**2 - 0.057638157095631 * (T[pos11] + 2) + 0.746900791092860)

    #Now make sure that rigidity is positive
    pos = np.nonzero(rigidity < 0)[0]
    if len(pos):
        rigidity[pos] = 1.0e6

    return rigidity
