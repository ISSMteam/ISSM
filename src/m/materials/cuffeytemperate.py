import numpy as np
import cuffey


def cuffeytemperate(temperature, waterfraction, stressexp):
    '''
    CUFFEYTEMPERATE - calculates ice rigidity as a function of temperature and waterfraction

   rigidity (in s^(1 / 3)Pa) is the flow law parameter in the flow law sigma = B * e(1 / 3)
   (Cuffey and Paterson, p75).
   temperature is in Kelvin degrees

   Usage:
      rigidity = cuffeytemperate(temperature, waterfraction, stressexp)
    '''

    if np.any(temperature < 0.):
        raise RuntimeError("input temperature should be in Kelvin (positive)")

    if (np.any(temperature.shape in waterfraction.shape)):
        error('input temperature and waterfraction should have same size!')

    if np.any(waterfraction < 0 | waterfraction > 1):
        error('input waterfraction should be between 0 and 1')

    rigidity = np.multiply(cuffey(temperature), (1 * np.ones(waterfraction.shape) + 181.25 * np.maximum(np.zeros(waterfraction.shape), np.minimum(0.01 * np.ones(waterfraction.shape), waterfraction)))**(-1 / stressexp))

    return rigidity
