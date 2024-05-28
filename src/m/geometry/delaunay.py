import numpy as np
from scipy.spatial import Delaunay # See also: https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.Delaunay.html?highlight=delaunay#scipy.spatial.Delaunay


def delaunay(*args):
    """DELAUNAY - Delaunay triangulation

    Replicates behaviour of MATLAB's 'delaunay'.

    Usage:
        # Creates a 2-D or 3-D Delaunay triangulation from the points in a 
        matrix P. The output DT is a three-column (for two dimensions) or 
        four-column (for three dimensions) matrix where each row contains the 
        row indices of the input points that make up a triangle or tetrahedron 
        in the triangulation.

        DT = delaunay(P)

        # Creates a 2-D Delaunay triangulation from the points in vectors x and 
        y.

        DT = delaunay(x, y)

        # Creates a 3-D Delaunay triangulation from the points in vectors x, y, 
        and z.

        DT = delaunay(x, y, z)

    Sources:
    - https://www.mathworks.com/help/matlab/ref/delaunay.html
    """

    nargs = len(args)

    if nargs == 1:
        points = args[0]
    elif nargs == 2:
        points = np.vstack((args[0], args[1])).T
    elif nargs == 3:
        # NOTE: Not tested, but it is assumed that 3-D triangulation would work.
        points = np.vstack((args[0], args[1], args[1])).T
    else:
        print(delaunay.__doc__)
        raise Exception('Wrong usage (see above)')

    simplices = Delaunay(points).simplices.astype(int) # NOTE: Need to covert to array of int so that it can be fetched properly by core

    return simplices
