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

def interpISMIP7AntarcticaOcn(*args): # {{{
    '''
    interpISMIP7AntarcticaOcn - interpolate chosen ISMIP7 ocean forcing to model

        Globus directory:
        AIS/
            CESM2-WACCM/
                historical/
                ssp126/
                ssp245/
                ssp585/
            obs/
            meltmip/
            SMBmip/
            parameterizations/

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
            md.basalforcings = interpISMIP7AntarcticaOcn(md,'cesm2-waccm','ssp126',[1995, 2100])
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
        datadir = '/data2/msmg/DATA/ISMIP7/AIS/'
    elif hostname == 'simba41':
        datadir = '/data04/Data/ISMIP7/AIS/'
    else:
        raise Exception('machine not supported yet, please provide your own path')

    # Search forcing files
    tf_file, so_file = search_forcing_file(datadir, modelname, scenario)

    # Define field name depending on modelname
    if modelname == 'obs':
        tf_name = 'tf'
        so_name = 'so'
    else:
        tf_name = 'tf'
        so_name = 'so'

    # Load TF and salinity data
    # First, load grid information. 
    nc = netCDF4.Dataset(tf_file[0],'r')
    x_n = nc['x'][:]
    y_n = nc['y'][:]
    # Python: dimension (time, z, y, x) for tf and so file
    z_data  = nc['z'][:]

    # Close...
    nc.close()

    tf_data   = []
    so_data   = []
    time_data = []
    print('   == loading Thermal Forcing (TF)')
    for i in range(len(tf_file)):
        nc_tf = netCDF4.Dataset(tf_file[i],'r')
        print(tf_file[i])
        tf_data = np.concatenate((tf_data, nc_tf[tf_name][:]),axis=0)
        try:
            time_data = np.concatenate((time_data, nc_tf['time']),axis=0)
        except:
            continue
        nc_tf.close()
    print('   == loading Salinity (SO)')
    for i in range(len(so_file)):
        print(so_file[i])
        nc_so = netCDF4.Dataset(so_file[i],'r')
        so_data = np.concatenate((so_data, nc_so[so_name][:]),axis=0)
        nc_so.close()

    # Corret time
    time_data = time_data/365 + 1850

    # Build tf and salinity array
    tf = []
    so = []
    if modelname:
        start_idx = 0
        final_idx = 1
        time = np.array([[1996]])
    elif modelname.lower() in ['cesm2-waccm']:
        start_idx = np.where(time_data == start_time)[0][0]
        final_idx = np.where(time_data == end_time)[0][0]
        time = np.arange(start_time, end_time+1)
    else:
        raise Exception('Error: Given' + modelname + ' is not supported.')

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
    cal_gamma, cal_delta_t = calibrated_parameters_ismip7()

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
    basalforcings.delta_t     = cal_delta_t
    basalforcings.tf_depths   = np.reshape(z_data,(1,-1))
    basalforcings.tf          = tf
    basalforcings.salinity    = so
    basalforcings.gamma       = cal_gamma

    print('Info: forcings cover ' + str(start_time) + ' to ' + str(end_time));

    return basalforcings
    # }}}

def search_forcing_file(datadir, modelname, scenario, start_time, end_time): # {{{
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

    start_time, end_time: int or float
        Start and fnal year for searching files

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

        tf_file = [tf_file]
        so_file = [so_file]
    elif modelname == 'cems2-waccm':
        tf_file_hist = glob.glob(os.path.join(dataset,'CESM2-WACCM','historical','ocean/tf/v3/tf*.nc'))
        tf_file_proj = glob.glob(os.path.join(dataset,'CESM2-WACCM',scenario,'ocean/tf/v3/tf*.nc'))

        tf_file_hist = np.sort(tf_file_hist)
        tf_file_proj = np.sort(tf_file_proj)

        so_file_hist = glob.glob(os.path.join(dataset,'CESM2-WACCM','historical','ocean/so/v3/so*.nc'))
        so_file_proj = glob.glob(os.path.join(dataset,'CESM2-WACCM',scenario,'ocean/so/v3/so*.nc'))

        so_file_hist = np.sort(so_file_hist)
        so_file_proj = np.sort(so_file_proj)

        # Merget file lists. Coerce list to array
        tf_file = np.array(tf_file_hist + tf_file_proj)
        so_file = np.array(so_file_hist + so_file_proj)

        # Choose specific year
        #NOTE:
        #File format: tf_AIS_CESM2-WACCM_historical_ocean_v3_2000-2009.nc
        years = np.arange(start_time, end_time+1)

        pos = np.zeros((len(tf_file),))
        for i in range(len(tf_file)):
            tmp_year = os.path.split(tf_file[i])[-1] # get file anme
            tmp_year = os.path.splitext(tmp_year)[0] # without extension
            tmp_year = tmp_year.split('_')[6] # select duration

            tmp_start = int(tmp_year.split('-')[0])
            tmp_end   = int(tmp_year.split('-')[1])

            # Now, check find file in years
            if np.any(years == tmp_start) | np.any(years == tmp_end):
                pos[i] = 1
        tf_file = tf_file[np.where(pos)[0])

        pos = np.zeros((len(so_file),))
        for i in range(len(so_file)):
            tmp_year = os.path.split(so_file[i])[-1] # get file anme
            tmp_year = os.path.splitext(tmp_year)[0] # without extension
            tmp_year = tmp_year.split('_')[6] # select duration

            tmp_start = int(tmp_year.split('-')[0])
            tmp_end   = int(tmp_year.split('-')[1])

            # Now, check find file in years
            if np.any(years == tmp_start) | np.any(years == tmp_end):
                pos[i] = 1
        so_file = so_file[np.where(pos)[0])
    else:
        raise Exception('Error: not implemented yet.')

    return tf_file, so_file
    # }}}

def calibrated_parameters_ismip7():# {{{
    """
    Explain
    -------
    Hard-coded optimized parameters for ismip7.
    
    Referneces
    ----------
    See notebook scripts at
    https://github.com/ismip/ismip7-antarctic-ocean-forcing/blob/main/parameterisations/parameter_selection_quadratic_example.ipynb
    """

    # Example
    #FIXME: unit for ISMIP7 protocol...
    yts = 31536000 # from md.constants.yts

    Kt = 7.5e-5*yts
    delta_t_basin = [-0.2,  -0.25, 0.15, 0.6 ,  0.1,
                    0.65, -0.2, -0.15, 0.8 ,  2.0,
                    0.55, -0.2,   0.5, 0.05, -0.2,
                    0.15]

    # Coerce list to array
    delta_t_basin = np.array(delta_t_basin)

    return Kt, deltaT_basin
    # }}}
