function basalforcings = interpISMIP7AntarcticaOcn(md, modelname, scenario, start_end)
	%interpISMIP6AntarcticaOcn - interpolate chosen ISMIP7 ocean forcing to model
	%
	%   Input:
	%     - md (model object)
	%     - modelname (string): name of the climate model and scenario
	%     - scenario  (string): scenario (i.e., ssp126, ssp585)
	%     - start_end (optional int array): two entry array of [start_year end_year]
	%
	%   Output:
	%     - basalforcings: prepared to be input directly into md.basalforcings
	%                      time series from 1995-2100
	%
	%   Examples:
	%      md.basalforcings = interpISMIP7AntarcticaOcn(md,'miroc-esm-chem_rcp8.5');
	%      md.basalforcings = interpISMIP7AntarcticaOcn(md,'miroc-esm-chem_rcp8.5', [2007 2050]);

	% Parse inputs
	if nargin==2 % for observation
		scenario   = '';
		start_time = 1996;
		end_time   = 1996;
	elseif nargin==3
		start_time = 1995;
		end_time = 2100;
	elseif nargin==4
		start_time = start_end(1);
		end_time = start_end(2);
	else
		error('no supported');
	end

	% Find appropriate directory
	% NOTE: data directory for ISMIP7 follows Globus repository...
	% Globus repository for ISMIP7.
	% https://app.globus.org/file-manager?origin_id=ccc9bbd2-4091-4e35-addd-eeb639cf5332&origin_path=%2FISMIP7%2F
	switch oshostname()
		case {'totten'}
			error('set default machine settings');
		case {'amundsen.thayer.dartmouth.edu'}
			error('set default machine settings');
		case {'simba00'}
			datadir='/data2/msmg/DATA/ISMIP7/AIS/';
		otherwise
			error('machine not supported yet, please provide your own path');
	end

	% Searching forcing files
	[tf_file, so_file] = search_forcing_file(datadir, modelname, scenario);

	%load TF data
	disp('   == loading TF');
	x_n     = double(ncread(tf_file,'x'));
	y_n     = double(ncread(tf_file,'y'));
	% dimension (x, y, z, time) for tf and so files.
	tf_data = double(ncread(tf_file,'tf')); 
	so_data = double(ncread(so_file,'tf')); % FIXME: really "tf" variable in "so" (salinity)?
	z_data  = double(ncread(tf_file,'z'));

	%Build tf cell array
	tf = cell(1,1,size(tf_data,3));
	so = cell(1,1,size(so_data,3));
	if modelname
		start_idx = 1;
		final_idx = 1;
		time = 1996; % set default starting time for observation.
	else
		start_idx = start_time - 1994;
		final_idx = end_time - 1994;
		time = start_time:end_time;
	end
	for i=1:size(tf_data,3)  %Iterate over depths
		disp(['   == Interpolating over depth ' num2str(i) '/' num2str(size(tf_data,3))]);
		
		temp_matrix_tf=[];
		temp_matrix_so=[];
		for ii=start_idx:final_idx %Iterate over time steps
			%temp_tfdata=InterpFromGridToMesh(x_n,y_n,tf_data(:,:,i,ii)',md.mesh.x,md.mesh.y,0);
			temp_data=InterpFromGrid(x_n,y_n,tf_data(:,:,i,ii)',md.mesh.x,md.mesh.y);
			temp_matrix_tf = [temp_matrix_tf temp_data];

			temp_data=InterpFromGrid(x_n,y_n,so_data(:,:,i,ii)',md.mesh.x,md.mesh.y);
			temp_matrix_so = [temp_matrix_so temp_data];
		end
		tf{:,:,i} = [temp_matrix_tf; time];
		so{:,:,i} = [temp_matrix_so; time];
	end

	clear temp_matrix_tf, temp_matrix_so;

	%load Delta and gamma data
	%deltatnc_median = fullfile(datadir,'parameterizations/coeff_gamma0_DeltaT_quadratic_non_local_median.nc');
	basin_datanc    = fullfile(datadir,'obs/ocean/IMBIE-basins/v3/IMBIE-basins_AIS_obs_ocean_v3.nc');
	%deltaT_median   = double(ncread(deltatnc_median,'deltaT_basin'));
	%gamma0_median   = double(ncread(deltatnc_median,'gamma0'));
	basinid_data    = double(ncread(basin_datanc,'basinNumber'));

	disp('   == Interpolating basin Id');
	num_basins = length(unique(basinid_data));
	%deltat_median = NaN(1,length(unique(basinid_data)));

	%for i=0:num_basins-1
	%	pos = find(basinid_data==i);
	%	deltat_temp = deltaT_median(pos);
	%	deltat_temp = deltat_temp(1);
	%	deltat_median(i+1) = deltat_temp;
	%end

	%Deal with basins ID
	x_el = mean(md.mesh.x(md.mesh.elements),2);
	y_el = mean(md.mesh.y(md.mesh.elements),2);
	basinid = InterpFromGrid(x_n,y_n,basinid_data',x_el, y_el, 'nearest')+1;

	%Set ISMIP7 basal melt rate parameters
	basalforcings            = basalforcingsismip7(md.basalforcings);
	basalforcings            = initialize(basalforcings,md);
	basalforcings.basin_id   = basinid;
	basalforcings.num_basins = num_basins;
	%basalforcings.delta_t    = deltat_median;
	basalforcings.tf_depths  = z_data';
	%basalforcings.gamma_0    = gamma0_median;
	basalforcings.tf         = tf;
	basalforcings.salinity   = so;

	disp(['Info: forcings cover ' num2str(start_time),' to ', num2str(end_time)]);
end

function [so_file, tf_file] = search_forcing_file(datadir, modelname, scenario)
	%{
	%Explain
	%	Return specific file names...
	%
	%Example
	%-------
	%.. code-block:: python
	%	fname = search_filenames
	%
	%Parameters
	%----------
	%datadir: str
	%
	%modelname: str
	%
	%scenario: str
	%
	%Returns
	%-------
	%tf_file, so_file: str
	%	thermal (tf_file) and salinity (so_file) forcing files, respectively.
	%}

	modelname = lower(modelname);
	switch modelname
		case 'obs'
			tf_file = fullfile(datadir,'obs/ocean/climatology/zhou_annual_06_nov/tf/v3/tf_AIS_obs_ocean_climatology_zhou_annual_06_nov_v3_1972-2024.nc');
			so_file = fullfile(datadir,'obs/ocean/climatology/zhou_annual_06_nov/so/v3/so_AIS_obs_ocean_climatology_zhou_annual_06_nov_v3_1972-2024.nc');
		case 'cesm2-waccm'
			tf_file = '';
			so_file = '';
		otherwise
			error('Error: not implemented yet.');
	end

	assert(exist(tf_file,'file')==2, ['Error: we cannot find filename: ' tf_file]);
	assert(exist(so_file,'file')==2, ['Error: we cannot find filename: ' so_file]);
end

function model_time_mapping(modelname, scenario, time_end)
	modelname = upper(modelname);
	switch modelname
		case 'CESM2-WACCM'
			historical=[[1850, 1859],...
							[1860, 1869],...
							[1870, 1879]];
		otherwise
			error('Error: not implemented yet.');
	end
end
