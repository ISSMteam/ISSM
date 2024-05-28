# module for inperpolating / smoothing data
import numpy as np
from scipy.interpolate import CloughTocher2DInterpolator, Rbf
from scipy.spatial import cKDTree
try:
    import matplotlib.pyplot as plt
except ImportError:
    print('could not import matplotlib, no plotting functions enabled. Set plotonly = False in function call')


def MeshSplineToMesh2d(x, y, data, xi, yi, tol=1e-6, fill_nans=False, **kwargs):  #{{{
    '''
    Piecewise cubic, C1 smooth, curvature-minimizing interpolant in 2D.
    The interpolant is guaranteed to be continuously differentiable,
    and the gradients are chosen such that the curvature of the interpolant
    is approximately minimized.

    Uses scipy.interpolate.CloughTocher2DInterpolator

    x, y:            data point coordinates
    data:            data to be interpolated (same length as x, y)
    xi, yi:        coordintes to interpolate data onto
    tol:            tolerance for gradient estimation (default 1e-6)
    fill_nans:    fill nan's (holes) in data using the spline fit?
    **kwargs:    optional keywork arguments:
                    maxiter: maximum iterations in gradient estimation

    Returns interpolated data at given x, y coordinates.

    Usage:
        interpdata = CloughToucher2d(x, y, data)

    Examples:
        interpdata = CloughToucher2d(md.mesh.x, md.mesh.y, data)
        interpdata = CloughToucher2d(md.mesh.x, md.mesh.y, data, tol = 1e-3, maxiter = 100)
    '''

    # unpack kwargs
    maxiter = kwargs.pop('maxiter', None)
    if 'maxiter' in kwargs:
        del kwargs['maxiter']
    if maxiter:
        assert type(maxiter) == int, 'error, maxiter should be an integer'
    assert len(kwargs) == 0, 'error, unexpected or misspelled kwargs'

    # create sub - vectors that just cover the limits of xi and yi
    # TODO x, y not necessarily a grid, so need a better definition of dx, dy (e.g. average element size)
    dx = 500
    dy = 500
    #dx = x[1] - x[0]
    #dy = y[1] - y[0]
    xlim = [min(xi) - dx, max(xi) + dx]
    ylim = [min(yi) - dy, max(yi) + dy]
    xflag = np.logical_and(x > xlim[0], x < xlim[1])
    yflag = np.logical_and(y > ylim[0], y < ylim[1])
    bothind = np.squeeze(np.where(np.logical_and(xflag, yflag))).astype(int)
    subdata = data[bothind]
    subx = x[bothind]
    suby = y[bothind]
    points = np.array([subx, suby]).T

    # mask out any nan's in the data and corresponding coordinate points
    mask = np.isnan(subdata)
    ind = np.nonzero(mask)[0]
    if len(ind) and fill_nans:
        print("        WARNING: filling nans using spline fit through good data points, which may or may not be appropriate. Check results carefully.")
    subdata = np.delete(subdata, ind)
    points = np.delete(points, ind, axis=0)

    if maxiter:
        spline = CloughTocher2DInterpolator(points, subdata, tol, maxiter=maxiter)
    else:
        spline = CloughTocher2DInterpolator(points, subdata, tol)

    interpdata = spline(xi, yi)

    if not fill_nans:
        # identify nan's in xi, yi using nearest neighbors
        xyinterp = np.dstack([xi, yi])[0]
        xg, yg = np.meshgrid(subx, suby)
        xydata = np.dstack([subx, suby])[0]
        tree = cKDTree(xydata)
        nearest = tree.query(xyinterp)[1]
        pos = np.nonzero(np.isnan(subdata[nearest]))
        interpdata[pos] = subdata[nearest][pos]

    return interpdata
    # }}}


def GridSplineToMesh2d(x, y, data, xi, yi, default_value=np.nan, plotonly=False, fill_nans=False):  #{{{
    '''
    python analog to InterpFromGridToMesh.  This routine uses
    scipy.interpolate.CloughTocher2dInterpolator to create a bivariate spline
    interpolation of the input data and then return values of the spline
    on the x, y coordinates of the model mesh.  The interpolant is piece-wise
    cubic, C1 smooth (continuously differentiable) and has approximately
    minimized curvature.  See "help(scipy.interpolate.CloughTocher2dInterpolator)"
    for more information on the routine.

    NOTE: this routine will not be appropriate if there are large holes (nan's) in
    the input data.  A non - spline interpolation scheme should be used in that case.

    x, y:                vectors defining the coordinates of the input data
    data:                2D array of input data
    xi, yi:            x and y coordinates to be interpolated onto
    default_value:    default value if points lie outside the convex hull of input
                        points (defaults to nan if not specified)
    plotonly:        plot the data to be interpolated using imshow (useful for
    fill_nans:        fill nan's (holes) in data using the spline fit?

    Usage:
        interpdata = GridToMesh(x, y, data, xi, yi, default_value = np.nan, plotonly = False, fill_nans = False)

    Examples:
        interpdata = GridToMesh(x_m, y_m, data, md.mesh.x, md.mesh.y, 0)
    '''

    if np.ndim(x) == 2:
        x = x.reshape(-1, )
    if np.ndim(y) == 2:
        y = y.reshape(-1, )
    if len(x) != data.shape[1] + 1 and len(x) != data.shape[1]:
        raise ValueError('x should have same length as ncols(data) or ncols(data) + 1')
    if len(y) != data.shape[0] + 1 and len(y) != data.shape[0]:
        raise ValueError('y should have same length as nrows(data) or nrows(data) + 1')

    # create sub - grid that just covers the limits of xi and yi
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    xlim = [min(xi) - dx, max(xi) + dx]
    ylim = [min(yi) - dy, max(yi) + dy]

    # TODO create grid differently depending on whether data is defined at x, y
    # or at the center of a grid cell with corner coordinates defined by xi, yi
    # create points array and flattened data array
    if len(x) == data.shape[1] and len(y) == data.shape[0]:
        print('        x, y taken to define the center of data grid cells')
        xind = np.nonzero(np.logical_and(x > xlim[0], x < xlim[1]))[0]
        yind = np.nonzero(np.logical_and(y > ylim[0], y < ylim[1]))[0]
        xg, yg = np.meshgrid(x[xind], y[yind])
        subdata = data[yind[0]:yind[-1] + 1, xind[0]:xind[-1] + 1]
    elif len(x) == data.shape[1] + 1 and len(y) == data.shape[0] + 1:
        print('        x, y taken to define the corners of data grid cells')
        xcenter = np.fromiter(((x[i] + x[i + 1]) / 2 for i in range(len(x) - 1)), np.float)
        ycenter = np.fromiter(((y[i] + y[i + 1]) / 2 for i in range(len(y) - 1)), np.float)
        xind = np.nonzero(np.logical_and(xcenter > xlim[0], xcenter < xlim[1]))[0]
        yind = np.nonzero(np.logical_and(ycenter > ylim[0], ycenter < ylim[1]))[0]
        xg, yg = np.meshgrid(xcenter[xind], ycenter[yind])
        subdata = data[yind[0]:yind[-1] + 1, xind[0]:xind[-1] + 1]
    else:
        raise ValueError('x and y have inconsistent sizes: both should have length ncols(data) / nrows(data) or ncols(data) + 1 / nrows(data) + 1')

    points = np.array([xg.ravel(), yg.ravel()]).T
    flatsubdata = subdata.ravel()

    if plotonly:
        plt.imshow(np.flipud(subdata), origin='upper')
        plt.show()
        return

    # mask out any nan's in the data and corresponding coordinate points
    mask = np.isnan(flatsubdata)
    ind = np.nonzero(mask)[0]
    if len(ind) and fill_nans:
        print("        WARNING: filling nans using spline fit through good data points, which may or may not be appropriate. Check results carefully.")
    goodsubdata = np.delete(flatsubdata, ind)
    goodpoints = np.delete(points, ind, axis=0)

    # create spline and index spline at mesh points
    spline = CloughTocher2DInterpolator(goodpoints, goodsubdata)
    interpdata = spline(xi, yi)

    if not fill_nans:
        # identify nan's in xi, yi using nearest neighbors
        xyinterp = np.dstack([xi, yi])[0]
        xydata = np.dstack([xg.ravel(), yg.ravel()])[0]
        tree = cKDTree(xydata)
        nearest = tree.query(xyinterp)[1]
        pos = np.nonzero(np.isnan(flatsubdata[nearest]))
        interpdata[pos] = flatsubdata[nearest][pos]

    return interpdata
    # }}}


def RadialInterp(x, y, data, xi, yi, **kwargs):  #{{{
    '''
    Interpolation using a radial basis function in 2 or 3 dimensions.
    Useful for smoothing input data after interpolation.

    Uses scipy.interpolate.Rbf

    x, y:            data point coordinates
    data:            data to be interpolated (same length as x, y)
    xi, yi:        coordinates to interpolate onto
    function:    form of radial basis function for interpolation:
                    'multiquadric': sqrt((r / self.epsilon)**2 + 1) (default)
                    'inverse': 1.0 / sqrt((r / self.epsilon)**2 + 1)
                    'gaussian': exp(-(r / self.epsilon)**2)
                    'linear': r
                    'cubic': r**3
                    'quintic': r**5
                    'thin_plate': r**2 * log(r)
    epsilon:        adjustable constant for scaling radial distance.  Defaults to
                    approximate average distance between nodes.
    smooth:        float > 0, adjusts the amount of smoothing applied.  Defaults to 0,
                    such that the function always passes through nodal points.
    z:                coordinate array if interpolating in 3 dimensions
    zi:            coordinate array if interpolating in 3 dimensions

    Usage:
        interpdata = RadialInterp(x, y, data,**kwargs)

    Examples:
        interpdata = RadialInterp(md.mesh.x, md.mesh.y, data)
        interpdata = RadialInterp(md.mesh.x, md.mesh.y, data, function = 'gaussian', epsilon = 100, smooth = 1)
    '''

    # unpack kwargs
    function = kwargs.pop('function', 'gaussian')
    if 'function' in kwargs:
        del kwargs['function']
    epsilon = kwargs.pop('epsilon', None)
    if 'epsilon' in kwargs:
        del kwargs['epsilon']
    smooth = kwargs.pop('smooth', 0)
    if 'smooth' in kwargs:
        del kwargs['smooth']
    z = kwargs.pop('z', None)
    if 'z' in kwargs:
        del kwargs['z']
    assert len(kwargs) == 0, 'error, unexpected or misspelled kwargs'

    if z:
        if epsilon:
            rbfi = Rbf(x, y, z, data, function=function, smooth=smooth, epsilon=epsilon)
        else:
            rbfi = Rbf(x, y, z, data, function=function, smooth=smooth)
        interpdata = rbfi(xi, yi, zi)
    else:
        if epsilon:
            rbfi = Rbf(x, y, data, function=function, smooth=smooth, epsilon=epsilon)
        else:
            rbfi = Rbf(x, y, data, function=function, smooth=smooth)
        interpdata = rbfi(xi, yi)

    return interpdata
    # }}}
