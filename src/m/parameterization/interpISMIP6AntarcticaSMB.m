function smb = interpISMIP6AntarcticaSMB(md,model_name)
%interpISMIP6AntarcticaSMB - interpolate chosen ISMIP6 atmospheric forcing to model
%
%   Input:
%     - md (model object)
%     - model_name (string): name of the climate model and scenario
%       - suppported options from /totten_1/ModelData/ISMIP6/Projections/AIS/Atmosphere_Forcing/
%             2.6 scenario             8.5 scenario
%             ---------------------------------------------
%             ccsm4_rcp2.6             ccsm4_rcp8.5
%                                      CESM2_ssp585 
%             CNRM_CM6_ssp126          CNRM_CM6_ssp585 
%                                      CNRM_ESM2_ssp585
%                                      CSIRO-Mk3-6-0_rcp85
%                                      HadGEM2-ES_rcp85
%             IPSL-CM5A-MR_rcp26       IPSL-CM5A-MR_rcp85
%             miroc-esm-chem_rcp2.6    miroc-esm-chem_rcp8.5
%             noresm1-m_rcp2.6         noresm1-m_rcp8.5
%
%   Output:
%     - smb: prepared to be input directly into md.smb
%            time series from 1995-2100
%
%   Examples:
%      md.smb = interpISMIP6AntarcticaSMB(md,'miroc-esm-chem_rcp8.5');

% Find appropriate directory
switch oshostname(),
	case {'totten'}
		path='/totten_1/ModelData/ISMIP6/Projections/AIS/Atmosphere_Forcing/';
	otherwise
		error('machine not supported yet, please provide your own path');
end

% search for smb files in the ISMIP climate model directory
rootname=[path model_name '/Regridded_2km/']; % root directory for the climate model files
smbclimfile=dir([rootname '*_2km_clim*_1995-2014.nc']); % climatology file if found
smbanomfile=dir([rootname '*_2km_anomaly*_1995-2100.nc']); % anomaly file if found

% throw error if files are not found, or if the file search is not unique
if length(smbclimfile)~=1 || length(smbanomfile)~=1
   error(['this path does not exist or is not unique under ' rootname]);
end

% save the full path of the found files
smbclimnc=[rootname smbclimfile.name];
smbanomnc=[rootname smbanomfile.name];

% load data from files
disp('   == loading TS and SMB climatology data');
lat                 = double(ncread(smbclimnc,'lat'));
lon                 = double(ncread(smbclimnc,'lon'));
smb_clim_data       = double(ncread(smbclimnc,'smb_clim'));
ts_clim_data        = double(ncread(smbclimnc,'ts_clim'));

disp('   == loading TS and SMB anomoly data');
smb_anomaly_data    = double(ncread(smbanomnc,'smb_anomaly'));
ts_anomaly_data     = double(ncread(smbanomnc,'ts_anomaly'));

%Create SMB and TS matrix
disp('   == Interpolating on model');
time = [1995:2100];
[x_n y_n]=ll2xy(lat(:,1),lon(:,1),-1);
y_n = x_n;
smb_clim = InterpFromGridToMesh(x_n,y_n,smb_clim_data',md.mesh.x,md.mesh.y,0);
ts_clim  = InterpFromGridToMesh(x_n,y_n,ts_clim_data',md.mesh.x,md.mesh.y,0);
temp_matrix_smb = []; temp_matrix_ts = [];
for i = 1:size(smb_anomaly_data,3)
	%SMB
	temp_smb        = InterpFromGridToMesh(x_n,y_n,smb_anomaly_data(:,:,i)',md.mesh.x,md.mesh.y,0);
	temp_smb        = temp_smb+smb_clim;
	temp_matrix_smb = [temp_matrix_smb temp_smb];
	%TS
	temp_ts         = InterpFromGridToMesh(x_n,y_n,ts_anomaly_data(:,:,i)',md.mesh.x,md.mesh.y,0);
	temp_ts         = temp_smb+smb_clim;
	temp_matrix_ts  = [temp_matrix_ts temp_ts];
	clear temp_smb; clear temp_ts;
end

% convert to m/yr
rhoi = md.materials.rho_ice;
temp_matrix_smb = temp_matrix_smb*(31556926/1000)*(1000/rhoi);

%Save Data (1995-2100)
smb = SMBforcing();
smb.mass_balance = [temp_matrix_smb ; time];

%What do we do with surface temp?
%md.miscellaneous.dummy.ts = [temp_matrix_ts ; time];
disp('Info: forcings cover 1995 to 2100');
end
