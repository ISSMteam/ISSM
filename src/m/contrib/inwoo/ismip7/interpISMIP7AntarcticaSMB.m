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
			path='/totten_1/ModelData/ISMIP6/Projections/AIS/Atmosphere_Forcing/';
		case {'amundsen.thayer.dartmouth.edu'}
			path='/local/ModelData/ISMIP6Data/Forcings2100/Atmosphere/';
		case {'simba00'}
			path='/data2/msmg/DATA/ISMIP7/AIS/';
		otherwise
			error('machine not supported yet, please provide your own path');
	end

	% Searching forcing files
	smb_file = search_forcing_file(datadir, modelname, scenario, start_time, end_time);

	% load data from files
	disp('   == loading SMB climatology data');
	x_n = double(ncread(smb_file{1},'x'));
	y_n = double(ncread(smb_file{1},'y'));	

	temp_matrix_smb = [];
	temp_matrix_time= [];
	for i = 1:length(smb_file)
		%NOTE: unit for acabf in netcdf file: kg m-2 s-1
		smb_data = double(ncread(smb_file{i},'acabf')); % dimension = (x,y,time)
		smb_data = smb_data/md.materials.rho_ice*md.constants.yts; % kg m-2 s-1 -> ice m yr-1

		temp_time = double(ncread(smb_file{i},'time'));
		temp_matrix_time = cat(1,temp_matrix_time, temp_time); % concatenate time series

		% Now, interpolate SMB 
		for j = 1:size(smb_data,3)
			temp_smb = InterpFromGridToMesh(x_n,y_n,smb_data(:,:,j)',md.mesh.x,md.mesh.y);

			% Concatenate dataset
			temp_matrix_smb = [temp_matrix_smb, temp_smb];
			clear temp_smb;
		end
	end

	clear smb_data, x_n, y_n;

	% convert days in year decimal
	%FIXME: standard calendar for time is 365 days in year (with noleap)?
	temp_matrix_time = temp_matrix_time/365 + 1850;

	% Save data
	smb = SMBforcing();
	smb.mass_balance = [temp_matrix_smb; temp_matrix_time];
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
			error('Not supported yet. We do not find any observation data set for SMB.';

		case 'cesm2-waccm'
			%FIXME: SDBN1 now is replaced with SDBN1-2000m or SDBN1-8000m. These search logic should be changed according to ISMIP7 repository.
			smb_file_hist = dir(fullfile(datadir,'CESM2-WACCM','historical','SDBN1/acabf/v2/acabf*.nc'));
			smb_file_proj = dir(fullfile(datadir,'CESM2-WACCM',scenario,'SDBN1/acabf/v2/acabf*.nc'));

			[~,pos]=sort({smb_file_hist.name});
			smb_file_hist = smb_file_hist(pos);
			[~,pos]=sort({smb_file_proj.name});
			smb_file_proj = smb_file_proj(pos);

			smb_file = cat(1,smb_file_hist, smb_file_proj);

			% Choose specific year
			%NOTE:
			%File format: acabf_AIS_CESM2-WACCM_historical_SDBN1_v2_2014.nc
			years = [start_time:1:end_time];

			pos = zeros(numel(smb_file),1);
			for i = 1:numel(smb_file)
				tmp_year = strsplit(smb_file(i).name,'_');
				tmp_year = tmp_year{7};
				tmp_year = strsplit(tmp_year,'.nc');
				tmp_year = int16(tmp_year{1});

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
