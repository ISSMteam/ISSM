import math

import numpy as np

from planetradius import *


def flagradiuselements(elements, x, y, z, lat0, long0, radius):  # {{{
    # get x0,y0,z0:
    R   = planetradius('earth')
    x0  = R * np.cos(np.deg2rad(lat0)) * np.cos(np.deg2rad(long0))
    y0  = R * np.cos(np.deg2rad(lat0)) * np.sin(np.deg2rad(long0))
    z0  = R * np.sin(np.deg2rad(lat0))

    distance    = np.array(math.sqrt(pow((x - x0), 2) + pow((y - y0), 2) + pow((z - z0), 2)))
    indices     = (distance <= radius * 1000).nonzero()

    #now that we know the indices, determine elements with own these indices:
    outelements = np.zeros(len(elements))
    for i in range(len(indices)):
        pos = np.array()
        pos = (elements == indices[i]).nonzero()
        for j in range(len(pos)):
            outelements[pos[j]] = 1
    outelements = outelements.nonzero()

    return outelements
# }}}
