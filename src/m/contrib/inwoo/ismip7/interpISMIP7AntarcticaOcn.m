function basalforcings = interpISMIP7AntarcticaOcn(md, modelname, scenario, start_end, varargin)
	%interpISMIP6AntarcticaOcn - interpolate chosen ISMIP7 ocean forcing to model
	%
	%   Globus directory:
	%   AIS/
	%      CESM2-WACCM/
	%         historical/
	%         ssp126/
	%         ssp245/
	%         ssp585/
	%	    obs/
	%      meltmip/
	%      SMBmip/
	%      parameterizations/
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
	%	Examples:
	%		# Get observation dataset
	%		md.basalforcings = interpISMIP7AntarcticaOcn(md,'obs')
	%
	%		md.basalforcings = interpISMIP7AntarcticaOcn(md,'cesm2-waccm','ssp126',[2007,2299]);
	%
	%  TODO:
	%  Do we really need to merge all forcings variables within single files using cdo? Or just search files in the data directory, which is synchronized with Globus?

	% Parse inputs
	if nargin==2 % for observation
		scenario   = 'obs';
		start_time = 2007;
		end_time   = 2007;
	elseif nargin==3 % for ESM model
		start_time = 1996;
		end_time   = 2299;
	elseif nargin>=4
		start_time = start_end(1);
		end_time = start_end(2);
	else
		error('not supported');
	end

	% Check consistency
	assert(any(strcmpi(scenario,{'obs','historical','ssp126','ssp585','ssp370'})));
	if any(strcmpi(scenario,{'historical'}))
		assert(end_time <= 2015,'Error: start_end(2) must be less than 2015 under ''historical'' scenario.');
	end
	if any(strcmpi(scenario,{'historical','ssp126','ssp585','ssp370'}))
		assert(end_time >= 2015,'Error: start_end(2) must be greater than 2015 under ''ssp'' scenarios.');
	end


	% Get options
	options = pairoptions(varargin{:});
	isdebug = getfieldvalue(options,'debug',false);

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
	[tf_file, so_file] = search_forcing_file(datadir, modelname, scenario, start_time, end_time);

	switch modelname
		case 'obs'
			%FIXME: Really salinity field name in observaiton is "tf"?
			tf_name = 'tf';
			so_name = 'so';
		otherwise
			tf_name = 'tf';
			so_name = 'so';
	end

	% Load TF and salinity data
	x_n     = double(ncread(tf_file{1},'x'));
	y_n     = double(ncread(tf_file{1},'y'));
	% dimension (x, y, z, time) for tf and so files.
	z_data  = double(ncread(tf_file{1},'z'));

	tf_data   = [];
	so_data   = [];
	time_data = [];
	disp('   == loading Thermal forcing (TF)');
	for i=1:numel(tf_file) 
		fprintf('      loading file %d/%d\r',i, numel(tf_file));
		tf_data = cat(4,tf_data,double(ncread(tf_file{i},tf_name)));
		try
			time_data = cat(1,time_data,double(ncread(tf_file{i},'time')));
		catch	
			continue
		end
	end
	fprintf('\n');

	disp('   == loading Salinity (SO)');
	for i=1:numel(so_file)
		fprintf('      loading file %d/%d\r',i, numel(tf_file));
		so_data = cat(4,so_data,double(ncread(so_file{i},so_name))); % FIXME: really "tf" variable in "so" (salinity)?
	end
	fprintf('\n');

	% Correct time
	time_data = time_data/365 + 1850;

	%Build tf and salinity cell array
	tf = cell(1,1,size(tf_data,3));
	so = cell(1,1,size(so_data,3));
	if strcmpi(modelname,'obs')
		start_idx = 1;
		final_idx = 1;
		time = 1996; % set default starting time for observation.
	elseif any(strcmpi(lower(modelname),{'cesm2-waccm'}))
		%Find start_idx and final_idx in given file
		start_idx = find(time_data == start_time);
		final_idx = find(time_data == end_time);
		if ~any(start_idx)
			start_idx = 1;
		end
		if ~any(final_idx)
			final_idx = size(time_data,1);
		end
		assert(length(start_idx) == 1);
		assert(length(final_idx) == 1);
		assert(start_idx < final_idx,	['Error: start_idx (=' num2str(start_idx) ') must be less than final_idx (=' num2str(final_idx) ').']);

		time = time_data(start_idx:final_idx);
		if size(time,1) ~= 1
			time = time';  % transpose to (1,ntime);
		end
	else
		error(['Error: Given ' modelanem ' is not supported.']);
	end

	%fprintf('   start_idx: %d, final_idx: %d\n',start_idx,final_idx);

	for i=1:size(tf_data,3)  %Iterate over depths
		disp(['   == Interpolating over depth ' num2str(i) '/' num2str(size(tf_data,3))]);
		
		temp_matrix_tf=[];
		temp_matrix_so=[];
		for nt=start_idx:final_idx %Iterate over time steps
			%temp_tfdata=InterpFromGridToMesh(x_n,y_n,tf_data(:,:,i,ii)',md.mesh.x,md.mesh.y,0);
			temp_data=InterpFromGrid(x_n,y_n,tf_data(:,:,i,nt)',md.mesh.x,md.mesh.y);
			temp_matrix_tf = [temp_matrix_tf, temp_data];

			temp_data=InterpFromGrid(x_n,y_n,so_data(:,:,i,nt)',md.mesh.x,md.mesh.y);
			temp_matrix_so = [temp_matrix_so, temp_data];
		end
		tf{:,:,i} = [temp_matrix_tf; time];
		so{:,:,i} = [temp_matrix_so; time];
	end

	% Clear memory: unused variables
	clear temp_matrix_tf, temp_matrix_so;
	clear so_data, tf_data;

	% TODO:
	% Wait calibrated dataset
	%load Delta and gamma data
	%deltatnc_median = fullfile(datadir,'parameterizations/coeff_gamma0_DeltaT_quadratic_non_local_median.nc');
	%deltaT_median   = double(ncread(deltatnc_median,'deltaT_basin'));
	%gamma0_median   = double(ncread(deltatnc_median,'gamma0'));
	[cal_gamma, cal_delta_t] = calibrated_parameters_ismip7();

	basin_datanc    = fullfile(datadir,'obs/ocean/IMBIE-basins/v3/IMBIE-basins_AIS_obs_ocean_v3.nc');
	basinid_data    = double(ncread(basin_datanc,'basinNumber'));

	disp('   == Interpolating basin Id');
	num_basins = length(unique(basinid_data));

	%Deal with basins ID
	x_el = mean(md.mesh.x(md.mesh.elements),2);
	y_el = mean(md.mesh.y(md.mesh.elements),2);
	basinid = InterpFromGrid(x_n,y_n,basinid_data',x_el, y_el, 'nearest')+1;

	%Set ISMIP7 basal melt rate parameters
	basalforcings            = basalforcingsismip7(md.basalforcings);
	basalforcings            = initialize(basalforcings,md);
	basalforcings.basin_id   = basinid;
	basalforcings.num_basins = num_basins;
	basalforcings.delta_t    = cal_delta_t;
	basalforcings.tf_depths  = z_data';
	basalforcings.tf         = tf;
	basalforcings.salinity   = so;
	basalforcings.gamma      = cal_gamma;

	disp(['Info: forcings cover ' num2str(start_time),' to ', num2str(end_time)]);
end

function [tf_file, so_file] = search_forcing_file(datadir, modelname, scenario, start_time, end_time) % {{{
	%{
	%Explain
	%-------
	%	Return specific file names...
	%
	%Example
	%-------
	%.. code-block:: matlab 
	%	[tf_file, so_file] = search_filenames(datadir, 'cesm2-waccm', 'ssp585')
	%
	%Parameters
	%----------
	%datadir: str
	%
	%modelname: str
	%
	%scenario: str
	%	Scenarios in CMIP6 (ssp126, ssp585, ssp370)
	%
	%start_time, end_time: int
	%	Start and final year for searching files
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
			so_file = fullfile(datadir,'obs/ocean/climatology/zhou_annual_06_nov/so/v4/so_AIS_obs_ocean_climatology_zhou_annual_06_nov_v4_1972-2024.nc');

			assert(isfile(tf_file),'Error: we cannot find tf_file for observation. %s',tf_file)
			assert(isfile(so_file),'Error: we cannot find so_file for observation. %s',so_file)

			tf_file = {tf_file};
			so_file = {so_file};
		case 'cesm2-waccm'
			tf_file_hist = dir(fullfile(datadir,'CESM2-WACCM','historical','ocean/tf/v3/tf*.nc'));
			tf_file_proj = dir(fullfile(datadir,'CESM2-WACCM',scenario,'ocean/tf/v3/tf*.nc'));

			[~,pos]=sort({tf_file_hist.name});
			tf_file_hist = tf_file_hist(pos);
			[~,pos]=sort({tf_file_proj.name});
			tf_file_proj = tf_file_proj(pos);

			so_file_hist = dir(fullfile(datadir,'CESM2-WACCM','historical','ocean/so/v3/so*.nc'));
			so_file_proj = dir(fullfile(datadir,'CESM2-WACCM',scenario,'ocean/so/v3/so*.nc'));

			[~,pos]=sort({so_file_hist.name});
			so_file_hist = so_file_hist(pos);
			[~,pos]=sort({so_file_proj.name});
			so_file_proj = so_file_proj(pos);

			if strcmpi(scenario,'historical')
				tf_file = tf_file_hist;
				so_file = so_file_hist;
			else % NOTE: Other scenarios.
				tf_file = cat(1,tf_file_hist, tf_file_proj);
				so_file = cat(1,so_file_hist, so_file_proj);
			end

			% Choose specific year
			%NOTE:
			%File format: tf_AIS_CESM2-WACCM_historical_ocean_v3_2000-2009.nc
			years = [start_time:1:end_time];

			pos = zeros(numel(tf_file),1);
			for i = 1:numel(tf_file)
				tmp_year = strsplit(tf_file(i).name,'_');
				tmp_year = tmp_year{7};
				tmp_year = strsplit(tmp_year,'.nc');
				tmp_year = tmp_year{1};

				tmp_start = strsplit(tmp_year,'-');
				tmp_start = str2num(tmp_start{1});
				tmp_end   = strsplit(tmp_year,'-');
				tmp_end   = str2num(tmp_end{2});

				% Now, check this find in years
				for tmp_year = tmp_start:tmp_end
					if any(tmp_year == years)
						pos(i) = 1;
						break
					end
				end
			end
			tf_file = tf_file(find(pos));

			pos = zeros(numel(so_file),1);
			for i = 1:numel(so_file)
				tmp_year = strsplit(so_file(i).name,'_');
				tmp_year = tmp_year{7};
				tmp_year = strsplit(tmp_year,'.nc');
				tmp_year = tmp_year{1};

				tmp_start = strsplit(tmp_year,'-');
				tmp_start = str2num(tmp_start{1});
				tmp_end   = strsplit(tmp_year,'-');
				tmp_end   = str2num(tmp_end{2});

				% Now, check this find in years
				for tmp_year = tmp_start:tmp_end
					if any(tmp_year == years)
						pos(i) = 1;
						break
					end
				end
			end
			so_file = so_file(find(pos));

			% Recover file names in cell.
			tf_file = fullfile({tf_file.folder},{tf_file.name});
			so_file = fullfile({so_file.folder},{so_file.name});

		otherwise
			error('Error: not implemented yet.');
	end

	assert(length(tf_file) > 0, ['Error: no thermal forcing files found.']);
	assert(length(so_file) > 0, ['Error: no salinity files found.']);
end % }}}

function [Kt,delta_t_basin]=calibrated_parameters_ismip7 % {{{
	%Explain
	%-------
	%Hard-coded optimized parameters for ismip7.
	%
	%Referneces
	%----------
	%See notebook scripts at
	%https://github.com/ismip/ismip7-antarctic-ocean-forcing/blob/main/parameterisations/parameter_selection_quadratic_example.ipynb

	% Example.
	%FIXME: unit for ISMIP7 protocol...
	yts = 31536000; % from md.constants.yts;

	iscase = 2;

	if iscase == 1
		Kt = 10^(-5.3737454)*yts;
		delta_t_basin = [0.77770618,  0.31441155,  0.80067001,  2.20641656,  1.098832,...
				1.51589658,  0.23850945,  0.79624342,  3.57994535,  4.86367174,...
				3.01984297,  2.89252678,  0.97019149,  1.29479108,  0.84029172,...
				1.48367211];
	elseif iscase == 2
		Kt = 10^(-4.799572805921373)*yts;
		delta_t_basin = [-0.9407537   0.58494678  0.45316308  0.08541776 -1.19998336 -0.66106441,...
			-0.98033621  0.27426227  0.88836344  1.83509204  0.58499364  0.00370923,...
			-0.88130295 -0.55646587  0.3053191   0.15185878];
	elseif iscase == 3
		Kt = 10^(-3.799572805921373)*yts;
		delta_t_basin = [-0.9407537   0.58494678  0.45316308  0.08541776 -1.19998336 -0.66106441,...
			-0.98033621  0.27426227  0.88836344  1.83509204  0.58499364  0.00370923,...
			-0.88130295 -0.55646587  0.3053191   0.15185878];
	else
		error('Error: not implemented yet.');
	end
end % }}}
