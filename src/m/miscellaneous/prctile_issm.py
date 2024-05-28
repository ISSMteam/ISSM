import numpy as np
from scipy.interpolate import interp1d


def prctile_issm(x, p, dim):
    # NumPy has no interpolation method that matches MATLAB's percentile function
    #y = np.percentile(x, p, dim, interpolation = 'higher')

    if len(np.shape(x)) > 2:
        raise RuntimeError('Number of dimensions {} not implemented.'.format(len(np.shape(x))))

    # presumably at least 1 input value has been given
    #    np.shape(integer) -> (), must be at least (1, )
    psize = np.shape(p) or (1, )
    if len(psize) > 1 and np.size(p, 1) > 1:
        p = p.T

    xsize = np.shape(x) or (1, )
    if dim == 2:
        x = x.T

    #  check for any NaN in any columns
    if not np.isnan(x).any():
        x = np.sort(x, axis=0)
        n = np.size(x, 0)

        #  branch based on number of elements
        if n > 1:
            #  set up percent values and interpolate
            xi = [((i + 0.5) * 100 / n) for i in range(n)]
            # scipy's interp1d returns a function
            y = interp1d(xi, x, axis=dim, bounds_error=False)
            y = y(p)

            #  fill in high and low values outside of interp range
            if p > xi[n - 1]:
                y = np.tile(x[n - 1, :], 1)
            if p < xi[0]:
                y = np.tile(x[0, :], 1)

        #  if one value, just copy it
        elif n == 1:
            if isinstance(p, int):
                y = x[0, :]
            else:
                y = np.tile(x[0, :], (len(p), 1))

        #  if no values, use NaN
        else:
            y = np.tile(float('NaN'), (np.size(p, 0), np.size(x, 0)))
    else:
        #  must loop over columns, since number of elements could be different
        y = np.zeros((np.size(p, 0), np.size(x, 1)))
        for j in range(np.size(x, 1)):
            #  remove any NaN and recursively call column
            y[:, j] = prctile_issm(x[np.where(not np.isnan(x[:, j]), j)], p)

    if (np.min(xsize) == 1 and len(xsize) > 1 and xsize[dim] > 1 and len(p) > 1 and psize[1] > 1) or (np.min(xsize) > 1 and dim == 2):
        y = y.T

    return y
