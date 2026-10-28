import datetime
import os

import numpy as np
from netCDF4 import Dataset
from InterpFromGridToMesh import InterpFromGridToMesh

def interpMIPkitGreenland(X, Y, string, ncfile=None, default_value=0): # {{{
    """INTERPMIPKITGREENLAND - interpolate ISMIP7 MIPkit Greenland datasets

    Available data:
        - mouginot_basins        (2004, doi:10.7280/D1WT11)
        - mask                   (2007: from BedMachine Greenland v6.6)
        - icemask_promice        (2022, Luetzenburg et al. 2025)
        - icemask_greene         (time series, nsidc-0793)
        - surface_gimp           (2007: from BedMachine Greenland v6.6)
        - surface_grimp          (2019, nsidc-0715)
        - thickness              (2007: from BedMachine Greenland v6.6)
        - bed                    (2007: from BedMachine Greenland v6.6)
        - geoid                  (EIGEN-6C4, Forste et al 2014)
        - vx_mosaic              (2010, nsidc-0670)
        - vy_mosaic              (2010, nsidc-0670)
        - vx_timeseries          (time series, nsidc-0478)
        - vy_timeseries          (time series, nsidc-0478)
        - dhdt_smith             (2019, Smith et al. 2020)
        - dhdt_khan              (time series, Khan et al. 2025)
        - geothermal_heat_flux1  (2021, Colgan et al. 2021)
        - geothermal_heat_flux2  (2018, Martos et al. 2018)

    For time series fields, the output has shape (len(X)+1, ntimes): the last
    row holds the time (decimal years). For other fields the output has the
    same shape as X.

    default_value is the value given to points outside of the data grid.

    Usage:
        dataout = interpMIPkitGreenland(X, Y, string, ncfile=None, default_value=0)

    Example:
        md.geometry.bed = interpMIPkitGreenland(md.mesh.x, md.mesh.y, 'bed')
        md.geometry.bed = interpMIPkitGreenland(md.mesh.x, md.mesh.y, 'bed', '../Data/GreenlandObsISMIP7-v1.3.nc')
    """

    verbose = True

    # Locate the netCDF file
    if ncfile is None:
        filename = 'GreenlandObsISMIP7-v1.3.nc'
        paths = [
            '/totten_1/ModelData/ISMIP7/MIPkit/' + filename,
            './Data/' + filename,
            './' + filename,
        ]
        for p in paths:
            if os.path.isfile(p):
                ncfile = p
                break
        else:
            raise FileNotFoundError('Could not find ' + filename + '. You can add the path to the list or provide its path as 4th argument')
    elif not os.path.isfile(ncfile):
        raise FileNotFoundError('File not found: ' + str(ncfile))

    X = np.asarray(X, dtype=np.float64)
    Y = np.asarray(Y, dtype=np.float64)

    # Field groups
    coarse_fields = {'dhdt_smith', 'geothermal_heat_flux1', 'geothermal_heat_flux2'}
    hires_fields = {'mouginot_basins', 'mask', 'icemask_promice',
                    'surface_gimp', 'surface_grimp', 'thickness', 'bed', 'geoid',
                    'vx_mosaic', 'vy_mosaic'}
    # field -> (x name, y name, time name)
    timeseries_fields = {
        'icemask_greene': ('x', 'y', 'greene_mask_time'),
        'dhdt_khan': ('x1km', 'y1km', 'khan_dhdt_time'),
        'vx_timeseries': ('x', 'y', 'vel_time'),
        'vy_timeseries': ('x', 'y', 'vel_time'),
    }

    with Dataset(ncfile, 'r') as nc:

        # Coordinates
        if verbose:
            print('   -- MIPkit: loading coordinates for ' + string)
        tdata = None
        if string in coarse_fields:      # non time series, 1km resolution
            xdata = readNCVariable(nc, 'x1km')
            ydata = readNCVariable(nc, 'y1km')
            istimeseries = False
        elif string in hires_fields:     # non time series, high resolution
            xdata = readNCVariable(nc, 'x')
            ydata = readNCVariable(nc, 'y')
            istimeseries = False
        elif string in timeseries_fields:
            xname, yname, tname = timeseries_fields[string]
            xdata = readNCVariable(nc, xname)
            ydata = readNCVariable(nc, yname)
            tdata = readNCVariable(nc, tname)
            istimeseries = True
        else:
            raise ValueError("data field '" + string + "' not supported")

        # Prepare subset based on provided X and Y (0-based inclusive indices)
        # x is assumed ascending, y descending (as in the MATLAB version)
        offset = 2
        nx, ny = xdata.size, ydata.size

        xmin, xmax = np.nanmin(X), np.nanmax(X)
        posx = np.flatnonzero(xdata <= xmax)
        last = posx[-1] if posx.size else nx - 1
        first = np.flatnonzero(xdata >= xmin)
        first = first[0] if first.size else nx - 1
        id1x = max(0, first - offset)
        id2x = min(nx - 1, last + offset)

        ymin, ymax = np.nanmin(Y), np.nanmax(Y)
        posy = np.flatnonzero(ydata >= ymin)
        last = posy[-1] if posy.size else ny - 1
        first = np.flatnonzero(ydata <= ymax)
        first = first[0] if first.size else ny - 1
        id1y = max(0, first - offset)
        id2y = min(ny - 1, last + offset)

        if verbose:
            print('   -- MIPkit: loading ' + string)
        # netCDF4 returns dims in file order: (y, x) or (time, y, x)
        if istimeseries:
            data = readNCVariable(nc, string, (slice(None), slice(id1y, id2y + 1), slice(id1x, id2x + 1)))
            data = np.moveaxis(data, 0, -1)          # -> (y, x, t)
        else:
            data = readNCVariable(nc, string, (slice(id1y, id2y + 1), slice(id1x, id2x + 1)))
        xdata = xdata[id1x:id2x + 1]
        ydata = ydata[id1y:id2y + 1]

    # InterpFromGridToMesh needs y in increasing order
    ydata = np.flipud(ydata)
    data = np.flipud(data)

    # Which interpolation type are we using? (as in the MATLAB version, this is only reported: InterpFromGridToMesh has no method argument)
    if string in {'mouginot_basins', 'mask', 'icemask_promice', 'icemask_greene'}:
        method = 'nearest'
    else:
        method = 'bilinear'

    # Go ahead and interpolate
    if verbose:
        print('   -- MIPkit: interpolating ' + string + ' (method: ' + method + ')')
    if istimeseries:
        if X.ndim > 1 and Y.ndim > 1 and min(X.shape) > 1:
            raise ValueError('X and Y should be vectors')
        tdata = date2decyear(tdata)
        dataout = np.zeros((X.size + 1, tdata.size))
        for i in range(tdata.size):
            slice_i = np.ascontiguousarray(data[:, :, i])
            dataout[:-1, i] = np.ravel(InterpFromGridToMesh(xdata, ydata, slice_i, X.ravel(), Y.ravel(), default_value))
            dataout[-1, i] = tdata[i]
    else:
        dataout = InterpFromGridToMesh(xdata, ydata, np.ascontiguousarray(data), X.ravel(), Y.ravel(), default_value)
        dataout = np.reshape(dataout, X.shape)

    return dataout
# }}}

def readNCVariable(nc, name, index=Ellipsis): # {{{
    """READNCVARIABLE - read a netCDF variable as float64 (fill values become NaN, like MATLAB ncread)

    Usage:
        data = readNCVariable(nc, name, index)
    """

    var = nc.variables[name]
    var.set_auto_maskandscale(True)
    out = np.ma.asarray(var[index]).astype(np.float64)
    return np.ma.filled(out, np.nan)
# }}}

def date2decyear(days): # {{{
    """DATE2DECYEAR - convert days since 1900-01-01 to decimal years

    Equivalent to ISSM's date2decyear.m (src/m/contrib/morlighem/ismip6) applied to
    datenum(datetime('1900-01-01') + days(days)): year + (days elapsed since Jan 1) / (days in that year)

    Usage:
        decyear = date2decyear(days)
    """

    base = datetime.datetime(1900, 1, 1)
    decyear = np.empty(len(days))
    for i, d in enumerate(days):
        dt = base + datetime.timedelta(days=float(d))
        start = datetime.datetime(dt.year, 1, 1)
        end = datetime.datetime(dt.year + 1, 1, 1)
        decyear[i] = dt.year + (dt - start) / (end - start)
    return decyear
# }}}