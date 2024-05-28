import numpy as np


def isnans(array):
    """
    ISNANS: figure out if an array is nan. wrapper to isnan from matlab which stupidly does not allow this test  for structures!

       Usage:    isnans(array)

          See also : ISNAN
    """

    if isinstance(array, (tuple, list, dict)):
        returnvalue = 0
    else:
        returnvalue = np.isnan(array)

    return returnvalue
