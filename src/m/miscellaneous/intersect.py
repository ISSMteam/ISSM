import numpy as np


def intersect(a, b):  # {{{
    """INTERSECT - Python implementation of MATLAB's 'intersect' function

    Usage:
        c = intersect(a, b) => Returns a list of values common to both 'a' and 'b', with no repetitions. 'c' is in sorted order.

        c, ia, ib = intersect(a, b) => Also returns index lists 'ia' and 'ib'.

    Sources:
    - https://www.mathworks.com/help/matlab/ref/double.intersect.html
    - https://stackoverflow.com/a/45645177
    """
    a_unique, ia = np.unique(a, return_index=True)
    b_unique, ib = np.unique(b, return_index=True)

    all_unique = np.concatenate((a_unique, b_unique))
    all_unique.sort()

    c = all_unique[:-1][all_unique[1:] == all_unique[:-1]]

    return c, ia[np.isin(a_unique, c)], ib[np.isin(b_unique, c)]
# }}}
