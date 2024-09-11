function frontalforcing = interpISMIP6GreenlandOcn(md,model_name)
%interpISMIP6GreenlandOcn - interpolate chosen ISMIP6 frontal forcing to model
%
%   This frontal melting implementation follows the parameterization proposed
%   by Slater et al. 2020 https://tc.copernicus.org/articles/14/985/2020/
%
%   Input:
%     - md (model object)
%     - model_name (string): name of the climate model and scenario
%     - suppported models:
%
%             2.6 scenario             8.5 scenario
%             ---------------------------------------------
%                                      access1-3_rcp8.5
%                                      csiro-mk3.6_rcp8.5
%                                      hadgem2-es_rcp8.5
%                                      ipsl-cm5-mr_rcp8.5
%             miroc-esm-chem_rcp2.6    miroc-esm-chem_rcp8.5
%                                      noresm1-m_rcp8.5      
%                                      ukesm1-cm6_ssp585
%   Output:
%     - undercutting: prepared to be input directly into md.basalforcings
%                      time series from 1950-2100
%
%   Examples:
%      md.frontalforcings = interpISMIP6GreenlandOcn(md,'miroc-esm-chem_rcp8.5');

% Find appropriate directory
switch oshostname(),
	case {'totten'}
		path='/totten_1/ModelData/ISMIP6/Projections/GrIS/Ocean_Forcing/Melt_Implementation/v4/';
	otherwise
		error('machine not supported yet, please provide your own path');
end

% search for thermal forcing file in the ISMIP climate model directory
rootname= [path '/' model_name ]; % root directory for the climate model files
tffile  = dir([rootname '/*oceanThermalForcing_v4.nc']); % thermal forcing file if found
qsgfile = dir([rootname '/*basinRunoff_v4.nc']); % subglacial discharge file if found

% throw error if file not found, or if the file search is not unique
if length(tffile)~=1
   error(['this path does not exist or is not unique under ' rootname]);
end
if length(qsgfile)~=1
	error(['this path does not exist or is not unique under ' rootname]);
end

% save the full path of the found files
tfnc  = [rootname '/' tffile.name];
qsgnc = [rootname '/' qsgfile.name];

%load TF data
disp('   == loading TF and Qsg');
x_n      = double(ncread(tfnc,'x'));
y_n      = double(ncread(tfnc,'y'));
tf_data  = double(ncread(tfnc,'thermal_forcing'));
time     = double(ncread(tfnc,'time'));
time_year= str2num(datestr(time+datenum(1900,1,1),'yyyy'));
qsg_data = double(ncread(qsgnc, 'basin_runoff'));

%Checks
for i=1:3
	assert(size(qsg_data,i)==size(tf_data,i));
end

%Interpolate on mesh
disp('   == Interpolating on model mesh');
TF_matrix  = zeros(md.mesh.numberofvertices, length(time_year));
Qsg_matrix = zeros(md.mesh.numberofvertices, length(time_year));
for i = 1:length(time_year)
	%TF
	TF_matrix(:,i)  = InterpFromGridToMesh(x_n,y_n, tf_data(:,:,i)',md.mesh.x,md.mesh.y, 0);
	%Qsg: change units from kg s-1 m-2 to m/day
	Qsg_matrix(:,i) = InterpFromGridToMesh(x_n,y_n, qsg_data(:,:,i)',md.mesh.x,md.mesh.y,0)*md.constants.yts/md.materials.rho_freshwater/365; 
end

%Remove negative values if need be
TF_matrix  = max(0, TF_matrix);
Qsg_matrix = max(0, Qsg_matrix);

%Set ISMIP6 frontal forcings
disp('   == Calculating undercutting melting rates');
A = 3e-4;
B = 0.15;
Alpha = 0.39;
Beta  = 1.18;
BED = -repmat(min(0,md.geometry.bed),[1 numel(time_year)]);

frontalforcing = frontalforcings();
frontalforcing.meltingrate = zeros(md.mesh.numberofvertices+1, numel(time_year));
frontalforcing.meltingrate(1:end-1,:) = (A*BED.*Qsg_matrix.^Alpha + B).*TF_matrix.^Beta *365;%Conversion from m/day to m/year
frontalforcing.meltingrate(end,:)     = time_year;

disp(['Info: forcings cover ' num2str(min(time_year)) ' to ' num2str(max(time_year))]);
