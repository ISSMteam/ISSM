import numpy as np


def nye(temperature, ice_type):
    """NYE - figure out the rigidity of ice (either CO2 or H2O) for a given 
    temperature rigidity (in s^(1/n)Pa) is the flow law parameter in the flow 
    law sigma=B*e(1/n) (Nye, p2000). Temperature is in Kelvin degrees.

    Usage:
        rigidity=nye(temperature,ice_type) % ice_type = 1: CO2 ice // ice_type = 2: H2O ice
    """

    # Declaring temperature and rigidity arrays
    if np.ndim(temperature) == 2:
        T = temperature.flatten()
    elif isinstance(temperature, float) or isinstance(temperature, int):
        T = np.array([temperature])
    else:
        T = temperature
    rigidity = np.zeros_like(T)

    # Beyond-melting-point cases
    if (ice_type == 1):
        for i in range(len(T)):
            if (200 < T[i] < 220):
                print('Warning: nye.py: CO2 ICE - POSSIBLE MELTING. Some temperature values are between 200K and 220K.')
            break
        if ((T >= 220).any()):
            print('Warning: nye.py: CO2 ICE - GUARANTEED MELTING. Some temperature values are beyond 220K.')
    elif (ice_type == 2) and ((T > 273.15).any()):
        print('Warning: nye.py: H2O ICE - GUARANTEED MELTING. Some temperature values are beyond 273.15K.')

    Rg = 8.3144598              # J mol^-1 K^-1

    if ice_type == 1:           # CO2 ice
        A_const = 1.0e13    # s^-1 MPa
        Q = 66900.              # J mol^-1
        n = 8.                  # Glen's exponent
    elif ice_type == 2:         # H2O ice
        A_const = 9e4       # s^-1 MPa
        Q = 60000.              #  J mol^-1
        n = 3.                  # Glen's exponent
    else:
        raise RuntimeError('Ice type not supported')

    # Arrhenius Law
    A = A_const * np.exp(-1 * Q / (T * Rg))  # s^-1 MPa
    rigidity = A**(-1 / n) * 1.0e6  # s^(1/n) Pa

    # Return output
    return rigidity
