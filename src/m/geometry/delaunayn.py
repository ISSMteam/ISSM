import numpy as np
from scipy.spatial import Delaunay # See also: https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.Delaunay.html?highlight=delaunay#scipy.spatial.Delaunay


def delaunayn(X, options=[]):
    """DELAUNAYN - N-D Delaunay triangulation

    Replicates behaviour of MATLAB's 'delaunayn'.

    Usage:
        # Computes a set of simplices such that no data points of X are 
        contained in any circumspheres of the simplices. The set of simplices 
        forms the Delaunay triangulation. X is an m-by-n array representing m 
        points in n-dimensional space. T is a numt-by-(n+1) array where each 
        row contains the indices into X of the vertices of the corresponding 
        simplex.

        T = delaunayn(X)

        # Specifies a cell array of options. The default options are:
        - {'Qt','Qbb','Qc'} for 2- and 3-dimensional input
        - {'Qt','Qbb','Qc','Qx'} for 4 and higher-dimensional input
        If options is [], the default options used. If options is {''}, no options are used, not even the default.

        T = delaunayn(X, options)

    Sources:
    - https://www.mathworks.com/help/matlab/ref/delaunayn.html
    - https://stackoverflow.com/questions/36604172/difference-between-matlab-delaunayn-and-scipy-delaunay
    """

    # Get dimensions of points
    ndims = X.shape[1]

    # Set up default options
    if options == []:
        options = 'Qt Qbb Qc' if ndims <= 3 else 'Qt Qbb Qc Qx'

    # Triangulate
    triangles = Delaunay(X, qhull_options=options).simplices

    # Strip the zero area/volume simplices that may have been created in the 
    # presence of degeneracy.
    keep = np.ones((len(triangles), ), dtype=bool)

    for i, t in enumerate(triangles):
        if abs(np.linalg.det(np.hstack((X[t], np.ones([1, ndims + 1]).T)))) < 1e-15:
            keep[i] = False

    triangles = triangles[keep].astype(int) # NOTE: Need to covert to array of int so that it can be fetched properly by core

    return triangles
