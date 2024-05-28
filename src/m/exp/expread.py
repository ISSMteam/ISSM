from collections import OrderedDict
import os.path

import numpy as np

import MatlabFuncs as m


def expread(filename):
    """expread - read a exp file and build a list of OrderedDicts

    This routine reads a file .exp and builds a list of OrderedDicts containing 
    the fields x and y corresponding to the coordinates, one for the filename 
    of the exp file, for the density, for the nodes, and a field closed to
    indicate if the domain is closed. The first argument is the .exp file to be 
    read and the second one (optional) indicate if the last point shall be read 
    (1 to read it, 0 not to).

    Usage:
        contours = expread(filename)

    Example:
        contours = expread('domainoutline.exp')
        contours = expread('domainoutline.exp')

    See Also:
    - expdoc
    - expwriteasvertices

    TODO:
    - Convert returned data structure from list of OrderedDict objects to list 
    of OrderedStruct objects (see also src/m/shp/shpread.py). Also, modify 
    handling of returned data structure in,
        - src/m/exp/expcoarsen.py
        - src/m/exp/expdisp.py
        - src/m/interp/SectionValues.py
        - src/m/mesh/bamg.py
    May also need to modify addressing in corresponding FetchData function, or
    create new one, in src/wrappers/ContoursToNodes/ContoursToNodes.cpp.
    """

    #some checks
    if not os.path.exists(filename):
        raise OSError("expread error message: file '%s' not found!" % filename)

    #initialize number of profile
    contours = []
    #open file
    fid = open(filename, 'r')
    #loop over the number of profiles
    while True:
        #update number of profiles
        contour = OrderedDict()
        #Get file name
        A = fid.readline()
        while A == '\n':
            A = fid.readline()
        if not A:
            break
        A = A.split(None, 1)
        if not (len(A) == 2 and m.strcmp(A[0], '##') and m.strncmp(A[1], 'Name:', 5)):
            break

        if len(A[1]) > 5:
            contour['name'] = A[1][5:-1]
        else:
            contour['name'] = ''

        #Get Icon
        A = fid.readline().split(None, 1)
        if not (len(A) == 2 and m.strcmp(A[0], '##') and m.strncmp(A[1], 'Icon:', 5)):
            break
        #Get Info
        A = fid.readline().split()
        if not (len(A) == 4 and m.strcmp(A[0], '#') and m.strcmp(A[1], 'Points')):
            break

        #Get number of nodes and density
        A = fid.readline().split()
        contour['nods'] = int(A[0])
        contour['density'] = float(A[1])

        #Get Info
        A = fid.readline().split()
        if not (len(A) == 5 and m.strcmp(A[0], '#') and m.strcmp(A[1], 'X') and m.strcmp(A[2], 'pos') and m.strcmp(A[3], 'Y') and m.strcmp(A[4], 'pos')):
            break
    #Get Coordinates
        contour['x'] = np.empty(contour['nods'])
        contour['y'] = np.empty(contour['nods'])
        for i in range(int(contour['nods'])):
            A = fid.readline().split()
            contour['x'][i] = float(A[0])
            contour['y'][i] = float(A[1])

    #Check if closed
        if (contour['nods'] > 1) and (contour['x'][-1] == contour['x'][0]) and (contour['y'][-1] == contour['y'][0]):
            contour['closed'] = True
        else:
            contour['closed'] = False

        contours.append(contour)
    #close file
    fid.close()
    return contours
