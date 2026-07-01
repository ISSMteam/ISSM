function smb = interpISMIP7AntarcticaSMB(md, modelname, scenario, start_end)
	%{
	%interpISMIP7AntarcticaSMB - interpolate chosen ISMIP7 SMB forcing to model.
	%
	%Example
	%-------
	%.. code-block:: matlab
	%	md.smb = interpISMIP7AntarcticaSMB(md,'CESM2-WACCM','ssp126',[1996, 2300]);
	%	md.smb = interpISMIP7AntarcticaSMB(md,'CESM2-WACCM','ssp585',[1996, 2300]);
	%
	%Inputs
	%------
	%md: ISSM model
	%
	%modelname: string
	%	CMIP6/CMIP7 model (i.e., CESM2-WACCM)
	%
	%scenario: string
	%	Scenario in CMIP6/7 model (i.e., ssp126, ssp585)
	%
	%start_end (optional int array): two entry array of [start_year end_year]
	%
	%Output
	%------
	%smb: SMBforcing class in ISSM.
	%}

	% Find appropriate directory
	switch oshostname(),
		case {'totten'}
			datadir='/totten_1/ModelData/ISMIP6/Projections/AIS/Atmosphere_Forcing/';
		case {'amundsen.thayer.dartmouth.edu'}
			datadir='/local/ModelData/ISMIP6Data/Forcings2100/Atmosphere/';
		case {'simba00'}
			datadir='/data2/msmg/DATA/ISMIP7/AIS/';
		otherwise
			error('machine not supported yet, please provide your own path');
	end

	if nargin == 3
		start_time = 1996;
		end_time   = 2300;
	elseif nargin == 4
		start_time = start_end(1);
		end_time   = start_end(2);
	end

	% Searching forcing files
	smb_file = search_forcing_file(datadir, modelname, scenario, start_time, end_time);

	% Load RACMO24p1_ERA5 dataset
	% Compute climatological mean value of SMB (Jan. 1995 - Dec. 2014 in Nowicki et al. (2020@TC))
	[smb_clim, ~] = interpRACMO24p1(md.mesh.x, md.mesh.y, 'smb','timebc',[1995,2014]);
	smb_clim = mean(smb_clim,2); % climatological mean value.

	% Load data from files
	disp('   == loading SMB anomaly data');
	x_n = double(ncread(smb_file{1},'x'));
	y_n = double(ncread(smb_file{1},'y'));	

	temp_matrix_smb_anon = [];
	temp_matrix_time= [];
	for i = 1:length(smb_file)
		fprintf('    processing file %d/%d \r',i,length(smb_file));

		%NOTE: unit for acabf in netcdf file: kg m-2 s-1
		smb_data = double(ncread(smb_file{i},'acabf-anomaly')); % dimension = (x,y,time)
		smb_data = smb_data/md.materials.rho_ice*md.constants.yts; % kg m-2 s-1 -> ice m yr-1

		%Load time data
		temp_time = double(ncread(smb_file{i},'time')); % time since year from current file...

		% configure out starting year of current file.
		temp_time_start = strsplit(smb_file{i},'_');
		temp_time_start = temp_time_start{7};
		temp_time_start = split(temp_time_start,'.nc');
		temp_time_start = str2int(temp_time_start{1});

		% convert days in year decimal
		%FIXME: standard calendar for time is 365 days in year (with noleap)?
		temp_time = temp_time/365 + temp_time_start;
		temp_matrix_time = cat(1,temp_matrix_time, temp_time); % concatenate time series

		% Now, interpolate SMB 
		for j = 1:size(smb_data,3)
			temp_smb_anon = InterpFromGridToMesh(x_n,y_n,smb_data(:,:,j)',md.mesh.x,md.mesh.y,NaN);

			% Concatenate dataset
			temp_matrix_smb_anon = [temp_matrix_smb_anon, temp_smb_anon];
			clear temp_smb_anon;
		end
	end


	clear smb_data;
	clear x_n, y_n;

	% Now, SMB = SMB_ref + SMB_anomaly
	temp_matrix_smb = repmat(smb_clim,1,size(temp_matrix_smb_anon,2)) + temp_matrix_smb_anon;	
	clear temp_matrix_smb_anon;

	% Check nan values
	if any(isnan(temp_matrix_smb))
		warning('There are NaN values in SMB forcing. Please check your data. We fill these NaN values with 0.0.');
		temp_matrix_smb(isnan(temp_matrix_smb)) = 0.0;
	end
	
	% Save data
	smb = SMBforcing();
	smb.mass_balance = [temp_matrix_smb; temp_matrix_time'];
end

function smb_file = search_forcing_file(datadir, modelname, scenario, start_time, end_time) % {{{
	%{
	%Explain
	%-------
	%	Return specific file names...
	%
	%Example
	%-------
	%.. code-block:: matlab 
	%	smb_file = search_filenames(datadir, 'cesm2-waccm', 'ssp585')
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
	%smb_file: cell with strings
	%	SMB forcing files, respectively.
	%}

	modelname = lower(modelname);

	switch modelname
		case 'obs'
			%FIXME: Assign smb file lists.
			error('Not supported yet. We do not find any observation data set for SMB.');

		case 'cesm2-waccm'
			%FIXME: SDBN1 now is replaced with SDBN1-2000m or SDBN1-8000m. These search logic should be changed according to ISMIP7 repository.
			smb_file_hist = dir(fullfile(datadir,'CESM2-WACCM','historical','SDBN1-8000m/acabf-anomaly/v2/acabf*.nc'));
			smb_file_proj = dir(fullfile(datadir,'CESM2-WACCM',scenario,'SDBN1-8000m/acabf-anomaly/v2/acabf*.nc'));

			[~,pos]=sort({smb_file_hist.name});
			smb_file_hist = smb_file_hist(pos);
			[~,pos]=sort({smb_file_proj.name});
			smb_file_proj = smb_file_proj(pos);

			if strcmpi(scenario,'historical')
				smb_file = smb_file_hist;
			else
				smb_file = cat(1,smb_file_hist, smb_file_proj);
			end

			% Choose specific year
			%NOTE:
			%File format: acabf_AIS_CESM2-WACCM_historical_SDBN1_v2_2014.nc
			years = [start_time:1:end_time];

			pos = zeros(numel(smb_file),1);
			for i = 1:numel(smb_file)
				tmp_year = strsplit(smb_file(i).name,'_');
				tmp_year = tmp_year{7};
				tmp_year = split(tmp_year,'.nc');
				tmp_year = str2int(tmp_year{1});

				% Now, check this find in years
				if any(years == tmp_year) 
					pos(i) = 1;
				end
			end
			smb_file = smb_file(find(pos));

			% Recover file names in cell.
			smb_file = fullfile({smb_file.folder},{smb_file.name});

		otherwise
			error('Error: not implemented yet.');
	end
end % }}}
