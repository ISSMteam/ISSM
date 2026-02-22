function basalforcings = interpISMIP6AntarcticaOcn(md, model_name, start_end)
%interpISMIP6AntarcticaOcn - interpolate chosen ISMIP6 atmospheric forcing to model
%
%   Input:
%     - md (model object)
%     - model_name (string): name of the climate model and scenario
%     - suppported models:
%
%             2.6 scenario             8.5 scenario
%             ---------------------------------------------
%                                      ccsm4_rcp8.5
%                                      cesm2_ssp585
%             cnrm-cm6-1_ssp126        cnrm-cm6-1_ssp585
%                                      cnrm-esm2-1_ssp585
%                                      csiro-mk3-6-0_rcp8.5
%                                      hadgem2-es_rcp8.5
%             ipsl-cm5a-mr_rcp2.6      ipsl-cm5a-mr_rcp8.5
%                                      miroc-esm-chem_rcp8.5
%             noresm1-m_rcp2.6         noresm1-m_rcp8.5
%                                      ukesm1-0-ll_ssp585
%     - start_end (optional int array): two entry array of [start_year end_year]
%
%   Output:
%     - basalforcings: prepared to be input directly into md.basalforcings
%                      time series from 1995-2100
%
%   Examples:
%      md.basalforcings = interpISMIP6AntarcticaOcn(md,'miroc-esm-chem_rcp8.5');
%      md.basalforcings = interpISMIP6AntarcticaOcn(md,'miroc-esm-chem_rcp8.5', [2007 2050]);

% Parse inputs
if nargin<3
   start_time = 1995;
   end_time = 2100;
elseif
   start_time = start_end(1);
   end_time = start_end(2);
else
   error('no supported');
end

% Find appropriate directory
switch oshostname(),
	case {'totten'}
	        path='/totten_1/ModelData/ISMIP6/Projections/AIS/Ocean_Forcing/';
        case {'amundsen.thayer.dartmouth.edu'}
                path='/local/ModelData/ISMIP6Data/Forcings2100/Ocean/';
	otherwise
		error('machine not supported yet, please provide your own path');
end

% search for thermal forcing file in the ISMIP climate model directory
rootname=[path model_name '/1995-2100/']; % root directory for the climate model files
tffile=dir([rootname '*_thermal_forcing_8km_x_60m.nc']); % thermal forcing file if found

% throw error if file not found, or if the file search is not unique
if length(tffile)~=1
   error(['this path does not exist or is not unique under ' rootname]);
end

% save the full path of the found files
tfnc=[rootname tffile.name];

%load TF data
disp('   == loading TF');
x_n     = double(ncread(tfnc,'x'));
y_n     = double(ncread(tfnc,'y'));
tf_data = double(ncread(tfnc,'thermal_forcing'));
z_data  = double(ncread(tfnc,'z'));

%Build tf cell array
time = start_time:end_time;
tf = cell(1,1,size(tf_data,3));
start_idx = start_time - 1994;
final_idx = end_time - 1994;
for i=1:size(tf_data,3)  %Iterate over depths
	disp(['   == Interpolating over depth ' num2str(i) '/' num2str(size(tf_data,3))]);
	
	temp_matrix=[];
	for ii=start_idx:final_idx %Iterate over time steps
		%temp_tfdata=InterpFromGridToMesh(x_n,y_n,tf_data(:,:,i,ii)',md.mesh.x,md.mesh.y,0);
		temp_tfdata=InterpFromGrid(x_n,y_n,tf_data(:,:,i,ii)',md.mesh.x,md.mesh.y);
		temp_matrix = [temp_matrix temp_tfdata];
	end
	tf{:,:,i} = [temp_matrix ; time];
end

%load Delta and gamma data
deltatnc_median = [path '/parameterizations/coeff_gamma0_DeltaT_quadratic_non_local_median.nc'];
basin_datanc    = [path '/imbie2/imbie2_basin_numbers_8km.nc'];
deltaT_median   = double(ncread(deltatnc_median,'deltaT_basin'));
gamma0_median   = double(ncread(deltatnc_median,'gamma0'));
basinid_data    = double(ncread(basin_datanc,'basinNumber'));

disp('   == Interpolating basin Id');
num_basins = length(unique(basinid_data));
deltat_median = NaN(1,length(unique(basinid_data)));

for i=0:num_basins-1
	pos = find(basinid_data==i);
	deltat_temp = deltaT_median(pos);
	deltat_temp = deltat_temp(1);
	deltat_median(i+1) = deltat_temp;
end

%Deal with basins ID
x_el = mean(md.mesh.x(md.mesh.elements),2);
y_el = mean(md.mesh.y(md.mesh.elements),2);
basinid = InterpFromGrid(x_n,y_n,basinid_data',x_el, y_el, 'nearest')+1;

%Set ISMIP6 basal melt rate parameters
basalforcings            = basalforcingsismip6(md.basalforcings);
basalforcings.basin_id   = basinid;
basalforcings.num_basins = num_basins;
basalforcings.delta_t    = deltat_median;
basalforcings.tf_depths  = z_data';
basalforcings.gamma_0    = gamma0_median;
basalforcings.tf         = tf;

disp(['Info: forcings cover ' num2str(start_time),' to ', num2str(end_time)]);
