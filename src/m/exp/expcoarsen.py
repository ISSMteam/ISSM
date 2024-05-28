import os.path
import numpy as np
from collections import OrderedDict
from expread import expread
from expwrite import expwrite


def expcoarsen(newfile, oldfile, resolution):
    """
    EXPCOARSEN - coarsen an exp contour

    This routine read an Argus file and remove points with respect to
    the resolution (in meters) given in input.

    Usage:
      expcoarsen(newfile, oldfile, resolution)

    Example:
       expcoarsen('DomainOutline.exp', 'Antarctica.exp', 4000)
    """

    #Some checks
    if not os.path.exists(oldfile):
        raise OSError("expcoarsen error message: file '%s' not found!" % oldfile)
    if os.path.exists(newfile):
        choice = eval(input('A file ' + newfile + ' already exists, do you want to modify it? (y / n)'))
        if choice not in 'y':
            print('no modification done ... exiting')
            return 0

    #Get exp oldfile
    contours = expread(oldfile)
    newcontours = []

    for contour in contours:
        numpoints = np.size(contour['x'])

        j = 0
        x = contour['x']
        y = contour['y']

        #stop if we have reached end of profile (always keep the last point)
        while j < numpoints - 1:

            #see whether we keep this point or not
            distance = np.sqrt((x[j] - x[j + 1])**2 + (y[j] - y[j + 1])**2)
            if distance < resolution and j < numpoints - 2:  #do not remove last point
                x = np.delete(x, j + 1, 0)
                y = np.delete(y, j + 1, 0)
                numpoints = numpoints - 1
            else:
                division = int(np.floor(distance / resolution) + 1)
                if division >= 2:
                    xi = np.linspace(x[j], x[j + 1], division)
                    yi = np.linspace(y[j], y[j + 1], division)

                    x = np.hstack((x[0:j + 1], xi[1:-1], x[j + 1:]))
                    y = np.hstack((y[0:j + 1], yi[1:-1], y[j + 1:]))

                    #update current point
                    j = j + 1 + division - 2
                    numpoints = numpoints + division - 2
                else:
                    #update current point
                    j = j + 1

        if np.size(x) > 1:
            #keep the (x, y) contour arond
            newcontour = OrderedDict()
            newcontour['nods'] = np.size(x)
            newcontour['density'] = contour['density']
            newcontour['x'] = x
            newcontour['y'] = y
            newcontours.append(newcontour)

    #write output
    expwrite(newcontours, newfile)
