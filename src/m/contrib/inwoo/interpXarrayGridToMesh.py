#!/usr/bin/env python3
import numpy as np
import xarray
import os, sys, platform

__all__ = ['interpXarrayGridToMesh']
def interpXarrayGridToMesh(fname, X, Y, varname=None, xname='lon', yname='lat', verbose:bool=False, method='linear', return_xarray:bool=True, selyear:list=None): # {{{
    '''interpXarray - interpolate netcdf file with xarray in python.

    Usage
    -----
    md = loadmodel('Models/Mesh.mat')
    fname = './pr_Amon_ACCESS-CM2_ssp126_r1i1p1f1_gn_201501-210012.nc'
    # ds - xarray output
    ds = interpXarrayGridToMesh(fname, md.mesh.long, md.mesh.lat, xname='lon', yname='lat')

    # Do interpolate Bedmachine dataset.
    md = loadmodel('Models/Mesh.mat')
    fname = 'BedMachineAntarctica-v3.nc'
    ds = interpXarrayGridToMesh(fname, md.mesh.x, md.mesh.y, xname='x',yname='y',varname='bed')

    Parameters
    ----------
    fname : str - netcdf file format.
    X,Y   : list or numpy.ndarray - do interpolation.

    method : str (default: linear) - interpolation method in xarray interpolation.
        Available methods
        linear, neareast, zero, sliner, quadratic, cubic, polynomial.
    xname : str (default: 'lon') - x coordinate name.
    yname : str (default: 'lat') - y coordinate name.
    verbose: bool (default: False) - show process.
    selyer: list (default: None) - select specific year range. i.g., use selyear = [1980, 2000] for selecting specific year between 1980, 2000.

    return_xarray: bool (default: True) - return results with xarray

    Returns
    -------
    output     - time series interpolate data.
    time       - return time.

    References
    ----------
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.interp.html
    '''

    # Check input format.
    if isinstance(fname,str):
        ds = xarray.load_dataset(fname)
    elif isinstance(fname,(xarray.core.dataset.Dataset, xarray.core.dataarray.DataArray)):
        ds = fname.copy()
    else:
        raise Exception('ERROR: Given data type is not available.')

    # Check interpolation method
    if not method in ['linear','nearest','zero','slinear','quadratic','cubic','polynomial']:
        raise Exception('ERROR: Given interpolation method (=%s) is not available. See the https://docs.xarray.dev/en/latest/generated/xarray.DataArray.interp.html for checking interpolation method'%(method))

    # Okay, check xname and yname in dimension of xarray file
    '''
    FutureWarning: The return type of `Dataset.dims` will be changed to return a set of dimension names in future, in order to be more consistent with `DataArray.dims`. To access a mapping from dimension names to lengths, please use `Dataset.sizes`.
    '''
    if isinstance(ds, (xarray.core.dataset.Dataset)):
        if xarray.__version__ > '2024.1.1':
            dims = list(ds.sizes.keys())
        else:
            dims = list(ds.dims.keys())
    elif isinstance(ds, (xarray.core.dataarray.DataArray)):
        if xarray.__version__ > '2024.1.1':
            dims = list(ds.sizes)
        else:
            dims = list(ds.dims)

    if not xname in dims:
        raise Exception(f'ERROR: Given xname (={xname}) is not defined in {dims}')
    elif not yname in dims:
        raise Exception(f'ERROR: Given yname (={yname}) is not defined in {dims}')

    # Initialize scatter array.
    xr = xarray.DataArray(X, dims='npts')
    yr = xarray.DataArray(Y, dims='npts')

    # First extract specific year
    if not (selyear is None):
        pos = [(yr in range(selyear[0], selyear[1]+1)) for yr in ds['time.year']]
        ds = ds.isel(time=pos)

    # Do interpolation
    if varname is not None:
        ds = ds[varname]
    interp = ds.interp(coords={xname: xr, yname:yr},
            method=method)

    # Return outputs
    if return_xarray:
        return interp
    else:
        return interp[varname], interp['time']
    # }}}

if __name__ == '__main__':
    import argparse

    # Do argument parsing.
    parser = argparse.ArgumentParser(prog='interpXarrayGridToMesh',
                                      description='Interpolate the netcdf file using xarray in python',
                                    )
    parser.add_argument('input',type=str,
                        help='Input file name. (e.g., input.nc)')
    parser.add_argument('variable',type=str,
                        help='Select specific variable in "input" file.')
    parser.add_argument('inputgrid',type=str,
                        help='Output grid format. (e.g.: input.xy)')
    parser.add_argument('output',type=str,
                        help='Output file. (e.g.: output.nc')
    parser.add_argument('-xname',type=str, default='lon',
                        help='name for x coordinate.')
    parser.add_argument('-yname',type=str, default='lat',
                        help='name for y coordinate.')
    parser.add_argument('-method',type=str, default='linear',
                        help='name for y coordinate.')
    parser.add_argument('-inputformat',type=str, default='txt',
                        )
    parser.add_argument('-v','--verbose',type=bool,default=False,
                        help='show process')
    args=parser.parse_args()

    # Okay, load x, y coordinates
    if args.inputformat == 'txt':
        xy = np.loadtxt(args.inputgrid)

        if args.verbose:
            print(xy)
    elif args.inputformat == 'mat':
        import mat73
        xy = mat73.loadmat(args.inputgrid)
        xy = np.vstack((xy['X'], xy['Y'])).T

    # Do interpolation
    ds = interpXarrayGridToMesh(args.input, xy[:,0], xy[:,1],
                                varname=args.variable,
                                xname=args.xname, yname=args.yname,
                                method=args.method,
                                verbose=args.verbose,
                                )

    if args.verbose:
        print('   interpXarrayGridToMesh: Export dataset.')
    ds.to_netcdf(args.output)
