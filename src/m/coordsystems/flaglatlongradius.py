import math

import numpy as np


def flaglatlogradius(lat, long, lat0, long0, radius):  # {{{
    '''
    FLAGLATLONGRADIUS - given a vector of lat, long, and a circle of radius 
                        degrees around lat0, long0, return the indices into 
                        lat, long that are within this circle.
                        lat and long should be between -90 and 90, and -180 and 
                        +180 respectively.
    '''

    #three cases, depending on whether our circle goes past the -180 +180 longitude line:
    if (long0 - radius) <= -180:
        for i in range(len(long)):
            if long[i] > 0:
                long[i] = long[i] - 360
    elif (long0 + radius) >= 180:
        for i in range(len(long)):
            if long[i] < 0:
                long[i] = long[i] + 360
    else:
        pass

    distance    = np.array(math.sqrt(pow((lat[i] - lat0[i]), 2) + pow((long[i] - long0[i]), 2)))
    indices     = (distance <= radius).nonzero() 

    return indices
# }}}
