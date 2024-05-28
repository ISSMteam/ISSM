import math

import numpy as np


def polyarea(x, y):  # {{{
    """POLYAREA - returns the area of the 2-D polygon defined by the vertices in 
    lists x and y
    
    Partial implementation of MATLAB's polyarea function. If x and y are 
    lists of the same length, then polyarea returns the scalar area of the 
    polygon defined by x and y.

    Usage:
        a = polyarea(x, y)

    Sources:
    - https://www.mathworks.com/help/matlab/ref/polyarea.html
    - https://stackoverflow.com/questions/24467972/calculate-area-of-polygon-given-x-y-coordinates
    - https://en.wikipedia.org/wiki/Shoelace_formula

    TODO:
    - Test that output falls within some tolerance of MATLAB's polyarea 
    function.
    """

    return 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))
# }}}
