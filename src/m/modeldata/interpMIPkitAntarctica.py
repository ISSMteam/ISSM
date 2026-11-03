import datetime
import os

import numpy as np
from netCDF4 import Dataset
from InterpFromGridToMesh import InterpFromGridToMesh

def interpMIPkitAntarctica(X, Y, string, ncfile=None, default_value=0): # {{{
    """INTERPMIPKITANTARCTICA - interpolate ISMIP7 MIPkit Antarctica datasets

    Available data:
        - mouginot_basins        (Rignot basins used by IMBIE-3, Rignot et al. 2011)
        - mask_bedmachine        (2015: BedMachine Antarctica v3, Morlighem et al. 2020)
        - surface_bedmachine     (2015: BedMachine Antarctica v3, Morlighem et al. 2020; REMA, Howat et al. 2019)
        - thickness_bedmachine   (2015: BedMachine Antarctica v3, Morlighem et al. 2020)
        - bed_bedmachine         (2015: BedMachine Antarctica v3, Morlighem et al. 2020)
        - mask_bedmap3           (2018: Bedmap3, Pritchard et al. 2025)
        - surface_bedmap3        (2015: Bedmap3, Pritchard et al. 2025; REMA, Howat et al. 2019)
        - thickness_bedmap3      (1966-2024: Bedmap3, Pritchard et al. 2025)
        - bed_bedmap3            (1966-2024: Bedmap3, Pritchard et al. 2025)
        - icemask_greene         (time series 1997-2021, Greene et al. 2022)
        - geoid                  (AntGG2021, Scheinert et al. 2024)
        - vx_mosaic              (1996-2016, nsidc-0484, Rignot et al. 2017)
        - vy_mosaic              (1996-2016, nsidc-0484, Rignot et al. 2017)
        - vx_timeseries          (time series, nsidc-0720, Mouginot et al. 2017)
        - vy_timeseries          (time series, nsidc-0720, Mouginot et al. 2017)
        - dhdt_smith             (2003-2019, Smith et al. 2020)
        - dhdt_cpom              (time series, ESA AIS CCI/CPOM, Shepherd et al. 2019)
        - geothermal_heat_flux1  (Stal et al. 2021, Aq1) 
        - geothermal_heat_flux2  (Losing & Ebbing 2021) 

    For time series fields, the output has shape (len(X)+1, ntimes): the last
    row holds the time (decimal years). For other fields the output has the
    same shape as X.

    default_value is the value given to points outside of the data grid.

    Usage:
        dataout = interpMIPkitAntarctica(X, Y, string, ncfile=None, default_value=0)

    Example:
        md.geometry.bed = interpMIPkitAntarctica(md.mesh.x, md.mesh.y, 'bed_bedmachine')
        md.geometry.bed = interpMIPkitAntarctica(md.mesh.x, md.mesh.y, 'bed_bedmachine', '../Data/AntarcticaObsISMIP7-v1.2.nc')
    """

    verbose = True

    # Locate the netCDF file
    if ncfile is None:
        filename = 'AntarcticaObsISMIP7-v1.2.nc'
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
    hires_fields = {'mouginot_basins',
                    'mask_bedmachine', 'surface_bedmachine', 'thickness_bedmachine', 'bed_bedmachine',
                    'mask_bedmap3', 'surface_bedmap3', 'thickness_bedmap3', 'bed_bedmap3',
                    'geoid', 'vx_mosaic', 'vy_mosaic'}
    # field -> (x name, y name, time name)
    timeseries_fields = {
        'icemask_greene': ('x', 'y', 'greene_mask_time'),
        'dhdt_cpom': ('x1km', 'y1km', 'cpom_dhdt_time'),
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
        offset = 2
        id1x, id2x = subsetIndices(xdata, np.nanmin(X), np.nanmax(X), offset)
        id1y, id2y = subsetIndices(ydata, np.nanmin(Y), np.nanmax(Y), offset)

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

    # InterpFromGridToMesh needs x and y in increasing order
    if ydata[0] > ydata[-1]:
        ydata = ydata[::-1]
        data = data[::-1, ...]
    if xdata[0] > xdata[-1]:
        xdata = xdata[::-1]
        data = data[:, ::-1, ...]

    # Which interpolation type are we using? (as for Greenland, this is only reported: InterpFromGridToMesh has no method argument)
    if string in {'mouginot_basins', 'mask_bedmachine', 'mask_bedmap3', 'icemask_greene'}:
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

def subsetIndices(coord, cmin, cmax, offset): # {{{
    """SUBSETINDICES - indices (0-based, inclusive) of a monotonic coordinate vector covering [cmin, cmax]

    The range is padded by offset points on each side. Works for ascending
    and descending coordinates.

    Usage:
        id1, id2 = subsetIndices(coord, cmin, cmax, offset)
    """

    n = coord.size
    descending = coord[0] > coord[-1]
    c = coord[::-1] if descending else coord

    below = np.flatnonzero(c <= cmax)
    last = below[-1] if below.size else 0
    above = np.flatnonzero(c >= cmin)
    first = above[0] if above.size else n - 1

    id1 = max(0, first - offset)
    id2 = min(n - 1, last + offset)
    if descending:
        id1, id2 = n - 1 - id2, n - 1 - id1
    return id1, id2
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