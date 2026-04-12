#!/usr/bin/env python3
import numpy as np
import os, sys, platform
import socket
try:
    import netCDF4
except:
    raise Exception('Error: netCDF4 is not installed. Install ''netCDF4'' module.')
try:
    import scipy
    import scipy.interpolate
except:
    raise Exception('Error: scipy is not installed. Install ''scipy'' module.')

from basalforcingsismip7 import basalforcingsismip7
from InterpFromGridToMesh import InterpFromGridToMesh

def interpISMIP7AntarcticaOcn(*args):
    '''
    interpISMIP7AntarcticaOcn - interpolate chosen ISMIP7 ocean forcing to model

        Input:
            - md (model object)
            - modelname (string): name of the climate model
            - scenario (string): scenario (e.g., ssp126, ssp370, ssp585)
            - start_end (optional int array): two array of [start_year, end_year]

        Output:
            - basalforcings: prepared to be input directly into md.basalforcings time series from 1995-2100.

        Example
            # Get observation dataset
            md.basalforcings = interpISMIP7AntarcticaOcn(md,'obs')

            # GCM forcings
            md.basalforcings = interpISMIP7AntarcticaOcn(md,'cesm2-waccm','ssp585',[1995, 2100])
    '''

    # Parse inputs
    if len(args) == 2:
        md, modelname = args
        scenario  = ''
        start_time= 1996
        end_time  = 1996
    elif len(args) == 3:
        md, modelname, scenario = args
        start_time = 1995
        end_time   = 2100
    elif len(args) == 4:
        md, modelname, scenario, start_end = args
    else:
        raise Exception('not supported.')

	# Find appropriate directory
	# NOTE: data directory for ISMIP7 follows Globus repository...
	# Globus repository for ISMIP7.
	# https://app.globus.org/file-manager?origin_id=ccc9bbd2-4091-4e35-addd-eeb639cf5332&origin_path=#2FISMIP7#2F
    hostname = socket.gethostname().lower().replace('-','')
    if hostname == 'totten':
        raise Exception('set default machine settings')
    elif hostname == 'amundsen.thayer.dartmouth.edu':
        raise Exception('set default machine settings')
    elif hostname == 'simba00':
        datadir='/data2/msmg/DATA/ISMIP7/AIS/'
    else:
        raise Exception('machine not supported yet, please provide your own path')

    # Search forcing files
    tf_file, so_file = search_forcing_file(datadir, modelname, scenario)

    # Load TF and salinity data
    nc_tf = netCDF4.Dataset(tf_file,'r')
    nc_so = netCDF4.Dataset(so_file,'r')

    x_n = nc_tf['x'][:]
    y_n = nc_tf['y'][:]
    # Python: dimension (time, z, y, x) for tf and so file
    tf_data = nc_tf['tf'][:]
    so_data = nc_so['tf'][:] # FIXME: really "tf" variable in "so" (salinity) ?
    z_data  = nc_tf['z'][:]
    if modelname == 'obs':
        #NOTE: observation dataset contains (z, y, x). observation dataset is required additional axis
        tf_data = tf_data[np.newaxis,:,:,:]
        so_data = so_data[np.newaxis,:,:,:]

    nc_tf.close()
    nc_so.close()
    del nc_tf, nc_so

    # Build tf and salinity array
    tf = []
    so = []
    if modelname:
        start_idx = 0
        final_idx = 1
        time = [[1996]]
    else:
        start_idx = start_time - 1994
        final_idx = end_time - 1994
        time = np.arange(start_time, end_time+1)

    for i in range(len(z_data)):
        print('   == Interpolating over depth ' + str(i+1) + '/' + str(len(z_data)))

        temp_matrix_tf=[]
        temp_matrix_so=[]
        for ii in range(start_idx,final_idx):
            temp_data=InterpFromGridToMesh(x_n,y_n,tf_data[ii,i,:,:],md.mesh.x,md.mesh.y,np.nan)
            temp_data=np.reshape(temp_data,(-1,1))
            temp_matrix_tf.append(temp_data[:])

            temp_data=InterpFromGridToMesh(x_n,y_n,so_data[ii,i,:,:],md.mesh.x,md.mesh.y,np.nan)
            temp_data=np.reshape(temp_data,(-1,1))
            temp_matrix_so.append(temp_data[:])

        temp_matrix_tf = np.concatenate(temp_matrix_tf,axis=1)
        temp_matrix_so = np.concatenate(temp_matrix_so,axis=1)

        tf.append(np.concatenate((temp_matrix_tf,time),axis=0))
        so.append(np.concatenate((temp_matrix_so,time),axis=0))

    del temp_matrix_tf, temp_matrix_so

    # TODO:
    # Wait calibrated dataset
    # load Delta and gamm data
    basin_datanc    = netCDF4.Dataset(os.path.join(datadir,'obs/ocean/IMBIE-basins/v3/IMBIE-basins_AIS_obs_ocean_v3.nc'),'r')
    basinid_data    = basin_datanc['basinNumber'][:]

    print('   == Interpolating basin Id')
    num_basins = len(np.unique(basinid_data))

    # Deal with basins ID
    x_el    = np.mean(md.mesh.x[md.mesh.elements-1],axis=1)
    y_el    = np.mean(md.mesh.y[md.mesh.elements-1],axis=1)
    # Use interpolator in "scipy.interpolate"
    interpolator = scipy.interpolate.RegularGridInterpolator((x_n,y_n),np.transpose(basinid_data),
                                                        method='nearest',
                                                        bounds_error=False,fill_value=np.nan)
    basinid = interpolator(np.vstack((x_el,y_el)).T)
    del interpolator

    # Set ISMIP7 basal melt rate parameters
    basalforcings             = basalforcingsismip7(md.basalforcings)
    basalforcings             = basalforcings.initialize(md)
    basalforcings.basin_id    = basinid
    basalforcings.num_basins  = num_basins
    basalforcings.tf_depths   = np.reshape(z_data,(1,-1))
    basalforcings.tf          = tf
    basalforcings.salinity    = so

    print('Info: forcings cover ' + str(start_time) + ' to ' + str(end_time));

    return basalforcings

def search_forcing_file(datadir, modelname, scenario):
    """
    Explain
    -------
        Return specific file names....

    Example
    -------
    .. code-block:: python
        tf_file, so_file = search_forcing_file(datadir, 'cesm2-waccm', 'ssp585')

    Parameters
    ----------
    datadir: str

    modelname: str

    scenario: str

    Returns
    -------
    tf_file, so_file: list
        thermal (tf_file) and salinity (so_file) forcing files, respectively.
    """

    assert(isinstance(modelname,str))
    modelname = modelname.lower()
    if modelname == 'obs':
        tf_file = os.path.join(datadir,'obs/ocean/climatology/zhou_annual_06_nov/tf/v3/tf_AIS_obs_ocean_climatology_zhou_annual_06_nov_v3_1972-2024.nc')
        so_file = os.path.join(datadir,'obs/ocean/climatology/zhou_annual_06_nov/so/v3/so_AIS_obs_ocean_climatology_zhou_annual_06_nov_v3_1972-2024.nc')
    elif modelname == 'cems2-waccm':
        raise Exception('Error: given %s is not supported yet.'%(modelname))
        tf_file = ''
        so_file = ''
    else:
        raise Exception('Error: not implemented yet.')

    assert(os.path.isfile(tf_file), 'Error: We cannot find filename: ' + tf_file)
    assert(os.path.isfile(so_file), 'Error: We cannot find filename: ' + so_file)

    return tf_file, so_file
