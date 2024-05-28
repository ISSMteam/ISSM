import numpy as np
from MatlabFuncs import strcmpi
from landcolor import landcolor
from seacolor import seacolor
from ibcao import ibcao


def demmap(ncolors, minZ, maxZ, colorscheme='dem'):
    '''DEMMAP - concatenate sea and land color deping on zmin and zmax

       Usage:
          cmap = demmap(n, zmin, zmax, colorscheme)

       Example:
          cmap = demmap(50, -300, 1200)
          cmap = demmap(50, -300, 1200, 'dem')
          cmap = demmap(50, -300, 1200, 'ibcao')
    '''

    if type(colorscheme) != str:
        raise RuntimeError('demmap: Error: optional argument "colorscheme" should be a string')

    # determine appropriate number of sea and land colors
    if minZ == maxZ:
        maxZ = minZ + 1

    cmn = minZ
    cmx = maxZ

    # determine appropriate number of sea and land colors
    if minZ >= 0:
        nsea = 0
        nland = ncolors
    elif maxZ <= 0:
        nland = 0
        nsea = ncolors
    else:
        # find optimal ratio of land to sea colors
        maxminratio = maxZ / abs(minZ)
        n1 = np.floor(ncolors / 2)
        n2 = np.ceil(ncolors / 2)
        if maxminratio > 1:
            sea = np.arange(1, n1 + 1)
            land = np.arange(ncolors - 1, n2 - 1, -1)
        else:
            land = np.arange(1, n1 + 1)
            sea = np.arange(ncolors - 1, n2 - 1, -1)

        ratio = land / sea
        errors = abs(ratio - maxminratio) / maxminratio
        indx = np.where(errors == min(errors))
        nsea = sea[indx]
        nland = land[indx]

    # determine color limits
        seaint = abs(minZ) / nsea
        landint = maxZ / nland
        if seaint >= landint:
            interval = seaint
        else:
            interval = landint

        cmn = -nsea * interval * (1 + 1e-9)  # zero values treated as land
        cmx = nland * interval

    if strcmpi(colorscheme, 'dem'):
        # concatenate and transpose to match matplotlib's colormap format
        cmap = np.concatenate((seacolor(nsea), landcolor(nland)**1.3), axis=1).T
    elif strcmpi(colorscheme, 'ibcao'):
        cmap = ibcao(nsea, nland)

    return cmap
