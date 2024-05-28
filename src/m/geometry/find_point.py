import numpy as np


def find_point(tabx, taby, pointx, pointy):
    """FIND_POINT - find closest point

    Find which point of the list (tabx, taby) is the closest to (pointx, pointy)

    Usage:
        f = find_point(tabx, taby, pointx, pointy)
    """

    # Compute distance between point and cloud of points
    distance = ((tabx - pointx) ** 2 + (taby - pointy) ** 2) ** 0.5 # NOTE: math.sqrt cannot be applied element-wise to a list/numpy.array

    # Find index of the minimum distance and return the first one only
    f = np.where(distance == np.min(distance))[0][0] # NOTE: numpy.where returns a tuple whose first element is a list, and we want the first element of that list

    return f
