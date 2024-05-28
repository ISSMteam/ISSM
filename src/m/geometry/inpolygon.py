import matplotlib.path as path
import numpy as np

def inpolygon(xq, yq, xv, yv):  # {{{
    """
    INPOLYGON - Returns points located inside polygonal region.

    Partial implementation of MATLAB's inpolygon function. Returns in 
    indicating if the query points specified by xq and yq are inside or on the 
    edge of the polygon area defined by xv and yv.

    Usage:
        in_polygon = inpolygon(xq, yq, xv, yv)

    Sources:
    - https://www.mathworks.com/help/matlab/ref/inpolygon.html
    - https://stackoverflow.com/questions/31542843/inpolygon-for-python-examples-of-matplotlib-path-path-contains-points-method/31543337#31543337
    """

    points = np.array((xq, yq)).T
    vertices = np.array((xv, yv)).T
    polygon = path.Path(vertices)
    in_polygon = polygon.contains_points(points)
    in_polygon = in_polygon.astype(int) # Convert from bool to int

    return in_polygon
# }}}
