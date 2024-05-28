%SNOWPACK class definition
%
%   Usage:
%      snowpack=snowpack();

classdef snowpack
	properties (SetAccess=public) 

		%first, the configuration fields, by category: 
		%snowpack:  %{{{
		snowpack_meas_tss = 0;
		snowpack_enforce_measured_snow_heights = 0;
		snowpack_sw_mode = 0;
		snowpack_incoming_longwave = 0;
		snowpack_height_of_wind_value = 0;
		snowpack_height_of_meteo_values = 0;
		snowpack_neutral = 0;
		snowpack_roughness_length = 0;
		snowpack_number_slopes = 0;
		snowpack_snow_redistribution = 0;
		snowpack_calculation_step_length = 0;
		snowpack_change_bc = 0;
		snowpack_thresh_change_bc = 0;
		snowpack_snp_soil = 0;
		snowpack_soil_flux = 0;
		snowpack_geo_heat = 0;
		snowpack_canopy = 0;
		%}}}
		%snowpackadvanced:  %{{{
		snowpackadvanced_variant = ''; % use 320 kg m-3 for fixed density
		snowpackadvanced_hn_density = '';
		%}}}
		%general:  %{{{
		general_pluginpath = '';
		general_buff_chunk_size = 0;
		general_buff_before = 0;
		%}}}
		%input {{{
		input_coordsys = '';
		input_coordparam = '';
		input_time_zone = 0;
		input_meteo = '';
		input_meteopath = '';
		input_station1 = '';
		input_snowfile1 = '';
		%}}}
		%output {{{
		output_coordsys = '';
		output_coordparam = '';
		output_time_zone = 0;
		output_meteopath = '';
		output_experiment = '';
		output_ts_write = 0;
		output_ts_start = 0;
		output_ts_days_between = 0;
		output_profile = '';
		output_prof_write = 0;
		output_prof_start = 0;
		output_prof_days_between = 0;
		%}}}
		%interpolations1d %{{{
		interpolations1d_window_size = 0; %that is 5 d and 2 h; 1 d = 86400
		interpolations1d_hnw_resample = '';
		interpolations1d_hs_resample = '';
		interpolations1d_tsg_resample = '';
		interpolations1d_rho_hn_resample = '';
		interpolations1d_vw_resample = '';
		interpolations1d_vw_args = '';
		%}}}
		%filters {{{
		filters={'TA::filter1',{'soft',[-20 10]}};
		filters=NaN;
		filter_values=NaN;

		filters_ta_filter1 = '';
		filters_ta_arg1 = NaN;
		filters_rh_filter1 = '';
		filters_rh_arg1 = NaN;
		filters_rh_filter2 = '';
		filters_rh_arg2 = NaN;
		filters_iswr_filter1 = '';
		filters_iswr_arg1 = NaN;
		filters_iswr_filter2 = '';
		filters_iswr_arg2 = NaN;
		filters_rswr_filter1 = '';
		filters_rswr_arg1 = NaN;
		filters_rswr_filter2 = '';
		filters_rswr_arg2 = NaN;

		%for ta between 190 and 280 k;
		filters_ilwr_filter1 = '';
		filters_ilwr_arg1 = NaN;
		filters_ilwr_filter2 = '';
		filters_ilwr_arg2 = NaN;
		filters_tss_filter1 = '';
		filters_tss_arg1 = NaN;
		filters_tsg_filter1 = '';
		filters_tsg_arg1 = NaN;
		filters_vw_filter1 = '';
		filters_vw_arg1 = NaN;
		filters_vw_filter2 = '';
		filters_vw_arg2 = NaN;
		%}}}

	end
	methods
		function self = snowpack(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					inputstruct=varargin{1};
					list1 = properties('snowpack');
					list2 = fieldnames(inputstruct);
					for i=1:length(list1)
						fieldname = list1{i};
						if ismember(fieldname,list2),
							self.(fieldname) = inputstruct.(fieldname);
						end
					end
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

		%snowpack:  %{{{
		self.snowpack_meas_tss = 1;
		self.snowpack_enforce_measured_snow_heights = 0;
		self.snowpack_sw_mode = 0;
		self.snowpack_incoming_longwave = 1;
		self.snowpack_height_of_wind_value = 12.;
		self.snowpack_height_of_meteo_values = 12.;
		self.snowpack_neutral = 0;
		self.snowpack_roughness_length = 0.002;
		self.snowpack_number_slopes = 1;
		self.snowpack_snow_redistribution = 1;
		self.snowpack_calculation_step_length = 15.0;
		self.snowpack_change_bc = 0;
		self.snowpack_thresh_change_bc = -1.0;
		self.snowpack_snp_soil = 0;
		self.snowpack_soil_flux = 0;
		self.snowpack_geo_heat = 0.06;
		self.snowpack_canopy = 0;
		%}}}
		%snowpackadvanced:  %{{{
		self.snowpackadvanced_variant = 'ANTARCTICA'; % use 320 kg m-3 for fixed density
		self.snowpackadvanced_hn_density = 'EVENT';
		%}}}
		%general:  %{{{
		self.general_pluginpath = '/usr/local/lib/meteoio/plugins/';
		self.general_buff_chunk_size = 90;
		self.general_buff_before = 1.5;
		%}}}
		%input {{{
		self.input_coordsys = 'ch1903';
		self.input_coordparam = 'null';
		self.input_time_zone = 8;
		self.input_meteo = 'smet';
		self.input_meteopath = './input';
		self.input_station1 = 'domec.smet';
		self.input_snowfile1 = 'domec.sno';
		%}}}
		%output {{{
		self.output_coordsys = 'ch1903';
		self.output_coordparam = 'null';
		self.output_time_zone = 8;
		self.output_meteopath = './output';
		self.output_experiment = 'smet';
		self.output_ts_write = 1;
		self.output_ts_start = 0.0;
		self.output_ts_days_between = 0.04166667;
		self.output_profile = 'ascii';
		self.output_prof_write = 1;
		self.output_prof_start = 0.0;
		self.output_prof_days_between = 0.04166667;
		%}}}
		%interpolations1d %{{{
		self.interpolations1d_window_size = 439200; %that is 5 d and 2 h; 1 d = 86400
		self.interpolations1d_hnw_resample = 'none';
		self.interpolations1d_hs_resample = 'linear';
		self.interpolations1d_tsg_resample = 'linear';
		self.interpolations1d_rho_hn_resample = 'none';
		self.interpolations1d_vw_resample = 'nearest_neighbour';
		self.interpolations1d_vw_args = 'extrapolate';
		%}}}
		%filters {{{
		self.filters_ta_filter1 = 'min_max';
		self.filters_ta_arg1 = [190 280];
		self.filters_rh_filter1 = 'min_max';
		self.filters_rh_arg1 = [0.01 1.2];
		self.filters_rh_filter2 = 'min_max';
		self.filters_rh_arg2 = {'soft' 0.01 1.0};
		self.filters_iswr_filter1 = 'min_max';
		self.filters_iswr_arg1 = [-10 1500];
		self.filters_iswr_filter2 = 'min_max';
		self.filters_iswr_arg2 = {'soft' 0 1500};
		self.filters_rswr_filter1 = 'min_max';
		self.filters_rswr_arg1 = [-10 1500];
		self.filters_rswr_filter2 = 'min_max';
		self.filters_rswr_arg2 = {'soft' 0 1500};

		%for ta between 190 and 280 k;
		self.filters_ilwr_filter1 = 'min_max';
		self.filters_ilwr_arg1 = [30 355];
		self.filters_ilwr_filter2 = 'min_max';
		self.filters_ilwr_arg2 = {'soft' 35 350};
		self.filters_tss_filter1 = 'min_max';
		self.filters_tss_arg1 = [180 275];
		self.filters_tsg_filter1 = 'min_max';
		self.filters_tsg_arg1 = [200 275];
		self.filters_vw_filter1 = 'min_max';
		self.filters_vw_arg1 = [-2 70];
		self.filters_vw_filter2 = 'min_max';
		self.filters_vw_arg2 = {'soft' 0 50};
		%}}}

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			%snowpack:  %{{{
			md=checkfield(md,'fieldname','snowpack.snowpack_meas_tss','values',[0 1]);
			md=checkfield(md,'fieldname','snowpack.snowpack_enforce_measured_snow_heights','values',[0 1]);
			md=checkfield(md,'fieldname','snowpack.snowpack_sw_mode','values',[0 1 2]);
			md=checkfield(md,'fieldname','snowpack.snowpack_incoming_longwave','values',[0 1]);
			md=checkfield(md,'fieldname','snowpack.snowpack_height_of_wind_value','>=',0);
			md=checkfield(md,'fieldname','snowpack.snowpack_height_of_meteo_values','>=',0);
			md=checkfield(md,'fieldname','snowpack.snowpack_neutral','values',[-1 0 1]);
			md=checkfield(md,'fieldname','snowpack.snowpack_roughness_length','>=',0);
			md=checkfield(md,'fieldname','snowpack.snowpack_number_slopes','values',[1 3 5 9]);
			md=checkfield(md,'fieldname','snowpack.snowpack_snow_redistribution','values',[0 1]);
			md=checkfield(md,'fieldname','snowpack.snowpack_calculation_step_length','>',0);
			md=checkfield(md,'fieldname','snowpack.snowpack_change_bc','values',[0 1]);
			md=checkfield(md,'fieldname','snowpack.snowpack_thresh_change_bc','<=',0);
			md=checkfield(md,'fieldname','snowpack.snowpack_snp_soil','values',[0 1]);
			md=checkfield(md,'fieldname','snowpack.snowpack_soil_flux','values',[0 1]);
			md=checkfield(md,'fieldname','snowpack.snowpack_geo_heat','>=',0);
			md=checkfield(md,'fieldname','snowpack.snowpack_canopy','values',[0 1]);
			%}}}
			%snowpackadvanced:  %{{{
			md=checkfield(md,'fieldname','snowpack.snowpackadvanced_variant','values',{'JAPAN','DEFAULT','ANTARCTICA'});
			md=checkfield(md,'fieldname','snowpack.snowpackadvanced_hn_density','values',{'PARAMETERIZED','EVENT','MEASURED'});
			%}}}
			%general:  %{{{
			md=checkfield(md,'fieldname','snowpack.general_buff_chunk_size','>',0);
			md=checkfield(md,'fieldname','snowpack.general_buff_before','>',0);
			%}}}
			%input {{{
			md=checkfield(md,'fieldname','snowpack.input_coordsys','values',{'CH1903','UTM','UPS','PROJ4','LOCAL'});
			md=checkfield(md,'fieldname','snowpack.input_coordparam','values','null');
			md=checkfield(md,'fieldname','snowpack.input_time_zone','>',-12,'<',12);
			md=checkfield(md,'fieldname','snowpack.input_meteo','values',{'BORMA','COSMO','GEOTOP','GRIB','GSN','IMIS','SMET','SNOWPACK'});
			md=checkfield(md,'fieldname','snowpack.input_meteopath','empty',1);
			md=checkfield(md,'fieldname','snowpack.input_station1 ','empty',1);
			md=checkfield(md,'fieldname','snowpack.input_snowfile1','empty',1);
			%}}}
			%output {{{
			md=checkfield(md,'fieldname','snowpack.output_coordsys','values',{'CH1903','UTM','UPS','PROJ4','LOCAL'});
			md=checkfield(md,'fieldname','snowpack.output_coordparam','values','null');
			md=checkfield(md,'fieldname','snowpack.output_time_zone','>',-12,'<',12);
			md=checkfield(md,'fieldname','snowpack.output_meteopath','empty',1);
			md=checkfield(md,'fieldname','snowpack.output_experiment','empty',1);
			md=checkfield(md,'fieldname','snowpack.output_ts_write','values',[0 1]);
			md=checkfield(md,'fieldname','snowpack.output_ts_start','>=',0);
			md=checkfield(md,'fieldname','snowpack.output_ts_days_between','>=',0);
			md=checkfield(md,'fieldname','snowpack.output_profile','values',{'ASCII','IMIS','ASCII IMIS'});
			md=checkfield(md,'fieldname','snowpack.output_prof_write','values',[0 1]);
			md=checkfield(md,'fieldname','snowpack.output_prof_start','>=',0);
			md=checkfield(md,'fieldname','snowpack.output_prof_days_between','>=',0);
			%}}}
			%interpolations1d %{{{
			md=checkfield(md,'fieldname','snowpack.interpolations1d_window_size','>',0);
			md=checkfield(md,'fieldname','snowpack.interpolations1d_hnw_resample','values',{'NONE','NEAREST_NEIGHBOUR','ACCUMULATE','LINEAR'});
			md=checkfield(md,'fieldname','snowpack.interpolations1d_hs_resample','values',{'NONE','NEAREST_NEIGHBOUR','ACCUMULATE','LINEAR'});
			md=checkfield(md,'fieldname','snowpack.interpolations1d_tsg_resample','values',{'NONE','NEAREST_NEIGHBOUR','ACCUMULATE','LINEAR'});
			md=checkfield(md,'fieldname','snowpack.interpolations1d_rho_hn_resample','values',{'NONE','NEAREST_NEIGHBOUR','ACCUMULATE','LINEAR'});
			md=checkfield(md,'fieldname','snowpack.interpolations1d_vw_resample','values',{'NONE','NEAREST_NEIGHBOUR','ACCUMULATE','LINEAR'});
			md=checkfield(md,'fieldname','snowpack.interpolations1d_vw_args','values',{'EXTRAPOLATE'});
			%}}}
			%filters {{{
			filter_values={'MIN_MAX','RATE_FILTER1','RATE_FILTER2','UNHEATED_RAIN_GAUGE_FILTER','WMO_UNDERCATCH_FILTER','WMO_UNDERCATCH_FILTER-SIMPLIFIED','UNVENTILLATED_TEMPERATURE_SENSOR','ADD_AN_OFFSET'};

			md=checkfield(md,'fieldname','snowpack.filters_ta_filter1','values',{filter_values});
			if strcmpi(md.snowpack.filters_ta_filter1,'MIN_MAX'), md=checkfield(md,'fieldname','snowpack.filters_ta_filter1','size',[1 NaN]); end
			md=checkfield(md,'fieldname','snowpack.filters_ta_arg1','values',{filter_values});
			if strcmpi(md.snowpack.filters_ta_arg1,'MIN_MAX'), md=checkfield(md,'fieldname','snowpack.filters_ta_arg1','size',[1 NaN]); end
			md=checkfield(md,'fieldname','snowpack.filters_rh_filter1','values',{filter_values});
			if strcmpi(md.snowpack.filters_rh_filter1,'MIN_MAX'), md=checkfield(md,'fieldname','snowpack.filters_rh_filter1','size',[1 NaN]); end
			md=checkfield(md,'fieldname','snowpack.filters_rh_arg1','values',{filter_values});
			if strcmpi(md.snowpack.filters_rh_arg1,'MIN_MAX'), md=checkfield(md,'fieldname','snowpack.filters_rh_arg1','size',[1 NaN]); end
			md=checkfield(md,'fieldname','snowpack.filters_rh_filter2','values',{filter_values});
			if strcmpi(md.snowpack.filters_rh_filter2,'MIN_MAX'), md=checkfield(md,'fieldname','snowpack.filters_rh_filter2','size',[1 NaN]); end
			md=checkfield(md,'fieldname','snowpack.filters_rh_arg2','values',{filter_values});
			if strcmpi(md.snowpack.filters_rh_arg2,'MIN_MAX'), md=checkfield(md,'fieldname','snowpack.filters_rh_arg2','size',[1 NaN]); end
			md=checkfield(md,'fieldname','snowpack.filters_iswr_filter1','values',{filter_values});
			if strcmpi(md.snowpack.filters_iswr_filter1,'MIN_MAX'), md=checkfield(md,'fieldname','snowpack.filters_iswr_filter1','size',[1 NaN]); end
			md=checkfield(md,'fieldname','snowpack.filters_iswr_arg1','values',{filter_values});
			if strcmpi(md.snowpack.filters_iswr_arg1,'MIN_MAX'), md=checkfield(md,'fieldname','snowpack.filters_iswr_arg1','size',[1 NaN]); end
			md=checkfield(md,'fieldname','snowpack.filters_iswr_filter2','values',{filter_values});
			if strcmpi(md.snowpack.filters_iswr_filter2,'MIN_MAX'), md=checkfield(md,'fieldname','snowpack.filters_iswr_filter2','size',[1 NaN]); end
			md=checkfield(md,'fieldname','snowpack.filters_iswr_arg2','values',{filter_values});
			if strcmpi(md.snowpack.filters_iswr_arg2,'MIN_MAX'), md=checkfield(md,'fieldname','snowpack.filters_iswr_arg2','size',[1 NaN]); end
			md=checkfield(md,'fieldname','snowpack.filters_rswr_filter1','values',{filter_values});
			if strcmpi(md.snowpack.filters_rswr_filter1,'MIN_MAX'), md=checkfield(md,'fieldname','snowpack.filters_rswr_filter1','size',[1 NaN]); end
			md=checkfield(md,'fieldname','snowpack.filters_rswr_arg1','values',{filter_values});
			if strcmpi(md.snowpack.filters_rswr_arg1,'MIN_MAX'), md=checkfield(md,'fieldname','snowpack.filters_rswr_arg1','size',[1 NaN]); end
			md=checkfield(md,'fieldname','snowpack.filters_rswr_filter2','values',{filter_values});
			if strcmpi(md.snowpack.filters_rswr_filter2,'MIN_MAX'), md=checkfield(md,'fieldname','snowpack.filters_rswr_filter2','size',[1 NaN]); end
			md=checkfield(md,'fieldname','snowpack.filters_rswr_arg2','values',{filter_values});
			if strcmpi(md.snowpack.filters_rswr_arg2,'MIN_MAX'), md=checkfield(md,'fieldname','snowpack.filters_rswr_arg2','size',[1 NaN]); end

			%for ta between 190 and 280 k;
			md=checkfield(md,'fieldname','snowpack.filters_ilwr_filter1','values',{filter_values});
			if strcmpi(md.snowpack.filters_ilwr_filter1,'MIN_MAX'), md=checkfield(md,'fieldname','snowpack.filters_ilwr_filter1','size',[1 NaN]); end
			md=checkfield(md,'fieldname','snowpack.filters_ilwr_arg1','values',{filter_values});
			if strcmpi(md.snowpack.filters_ilwr_arg1,'MIN_MAX'), md=checkfield(md,'fieldname','snowpack.filters_ilwr_arg1','size',[1 NaN]); end
			md=checkfield(md,'fieldname','snowpack.filters_ilwr_filter2','values',{filter_values});
			if strcmpi(md.snowpack.filters_ilwr_filter2,'MIN_MAX'), md=checkfield(md,'fieldname','snowpack.filters_ilwr_filter2','size',[1 NaN]); end
			md=checkfield(md,'fieldname','snowpack.filters_ilwr_arg2','values',{filter_values});
			if strcmpi(md.snowpack.filters_ilwr_arg2,'MIN_MAX'), md=checkfield(md,'fieldname','snowpack.filters_ilwr_arg2','size',[1 NaN]); end
			md=checkfield(md,'fieldname','snowpack.filters_tss_filter1','values',{filter_values});
			if strcmpi(md.snowpack.filters_tss_filter1,'MIN_MAX'), md=checkfield(md,'fieldname','snowpack.filters_tss_filter1','size',[1 NaN]); end
			md=checkfield(md,'fieldname','snowpack.filters_tss_arg1','values',{filter_values});
			if strcmpi(md.snowpack.filters_tss_arg1,'MIN_MAX'), md=checkfield(md,'fieldname','snowpack.filters_tss_arg1','size',[1 NaN]); end
			md=checkfield(md,'fieldname','snowpack.filters_tsg_filter1','values',{filter_values});
			if strcmpi(md.snowpack.filters_tsg_filter1,'MIN_MAX'), md=checkfield(md,'fieldname','snowpack.filters_tsg_filter1','size',[1 NaN]); end
			md=checkfield(md,'fieldname','snowpack.filters_tsg_arg1','values',{filter_values});
			if strcmpi(md.snowpack.filters_tsg_arg1,'MIN_MAX'), md=checkfield(md,'fieldname','snowpack.filters_tsg_arg1','size',[1 NaN]); end
			md=checkfield(md,'fieldname','snowpack.filters_vw_filter1','values',{filter_values});
			if strcmpi(md.snowpack.filters_vw_filter1,'MIN_MAX'), md=checkfield(md,'fieldname','snowpack.filters_vw_filter1','size',[1 NaN]); end
			md=checkfield(md,'fieldname','snowpack.filters_vw_arg1','values',{filter_values});
			if strcmpi(md.snowpack.filters_vw_arg1,'MIN_MAX'), md=checkfield(md,'fieldname','snowpack.filters_vw_arg1','size',[1 NaN]); end
			md=checkfield(md,'fieldname','snowpack.filters_vw_filter2','values',{filter_values});
			if strcmpi(md.snowpack.filters_vw_filter2,'MIN_MAX'), md=checkfield(md,'fieldname','snowpack.filters_vw_filter2','size',[1 NaN]); end
			md=checkfield(md,'fieldname','snowpack.filters_vw_arg2','values',{filter_values});
			if strcmpi(md.snowpack.filters_vw_arg2,'MIN_MAX'), md=checkfield(md,'fieldname','snowpack.filters_vw_arg2','size',[1 NaN]); end

			%}}}
		end % }}}
		function disp(self) % {{{

			disp(sprintf('   Snowpack solution parameters:'));
			disp(sprintf('\n	%s','Snowpack parameters:')); % {{{
			fielddisplay(self,'snowpack_meas_tss',{'A measured surface temperature is available and can be reliably ','used for various consistency tests (it needs to be set to true if enabling CHANGE_BC) (0 or 1)'});
			fielddisplay(self,'snowpack_enforce_measured_snow_heights','Input mode by which a measurement of snow depth is used to drive the snow cover mass balance (0 or 1)');
			fielddisplay(self,'snowpack_sw_mode',{'Define the shortwave radiation input:',...
				'0 Incoming shortwave radiation is measured and albedo estimated by the model',...
				'1 Reflected shortwave radiation is available as input and albedo is estimated by the model (IMIS standard)',...
				'2 Incoming and reflected shortwave radiation are both measured and the albedo is estimated from both measurements subject to plausibility checks.'});
			fielddisplay(self,'snowpack_incoming_longwave','Use the provided incoming long wave on the virtual slopes? (0 or 1)');
			fielddisplay(self,'snowpack_height_of_wind_value',{'The instrument height (or model layer height) for wind input data; note that height ',...
				'is above ground for a standard SNOWPACK application but above surface (snow or ground) for Alpine3D applications '});
			fielddisplay(self,'snowpack_height_of_meteo_values',{'The instrument height (or model layer height) for meteorological input data except for wind,',...
				'which may be at a different height; note that height is above ground for a standard SNOWPACK ',...
				'application but above surface (snow or ground) for Alpine3D applications. '});
			fielddisplay(self,'snowpack_neutral',{'Select the atmospheric stability correction model:',...
				'-1 use a simplified Richardson number stability correction',...
				'0 assume standard Monin-Obukhov bulk formulation for surface exchange iteration with Paulson, Stearns and Weidner (can be used with BC_CHANGE=0)',...
				'1 force Monin-Obukhov formulation to assume neutral conditions regardless of the actual stratification; it has been shown to work well in ',...
				'complex terrain settings. It should be used with BC_CHANGE=1, i.e., Dirichlet /* but also is recommended with Neumann b.c., i.e., BC_CHANGE=0.'});
			fielddisplay(self,'snowpack_roughness_length',{'Aerodynamic roughness length as a parameter for the Monin-Obukhov bulk formulation;',...
				'A typical value for complex terrain is 0.01 m and for snow covered flat sites 0.001 m. '});
			fielddisplay(self,'snowpack_number_slopes',{'Based on meteorological input from a (flat field) automatic station or numerical weather model,',...
				'up to 8 expositions can be calculated in addition to the flat field if the corresponding *.sno files are provided. For example,',...
				'if you provide a flat field *.snow file (mandatory), which is named KLO3.sno and you want 4 slopes to be calculated the corresponding',...
				'slope files should be named KLO21.sno, ...,KLO24.sno '});
			fielddisplay(self,'snowpack_snow_redistribution',{'Specifies if redistribution of snow is allowed from (upwind) expositions to lee slopes.',...
				'In case just the flat field is calculated, snow erosion is enabled but only for "ENFORCE_MEASURED_SNOW_HEIGHTS".'});
				fielddisplay(self,'snowpack_calculation_step_length',{'Internal time step (in minutes) used for model simulation. Please note that this MUST ',...
				'be the same as HNW::accumulate (the latter being in seconds) if re-acumulating precipitation, otherwise it would lead to wrong results.'});
			fielddisplay(self,'snowpack_change_bc',{'Use measured surface temperature as Dirichlet temperature BC for sub-freezing snowpacks and switch to ',...
			'Neumann only for melting snowpacks. If set to false, assumes Neumann boundary conditions.'});
			fielddisplay(self,'snowpack_thresh_change_bc','Threshold value (small number below freezing), which switches from Dirichlet to Neumann BCs if CHANGE_BC is selected');
			fielddisplay(self,'snowpack_snp_soil','Soil layers as defined by the *.sno files are included in the simulation');
			fielddisplay(self,'snowpack_soil_flux','Assume that the lower temperature boundary condition is given by GEO_HEAT (Neumann) and not by a measured temperature');
			fielddisplay(self,'snowpack_geo_heat','Constant geothermal heat flux at great) depth W m-2): Lower flux boundary condition for temperature equation if BC is Neumann');
			fielddisplay(self,'snowpack_canopy','Switch to tell the model that canopy is present (note that Canopy parameters should then be provided in the *.sno file)');
			% }}}
			disp(sprintf('\n	%s','Snowpackadvanced parameters:')); % {{{
			fielddisplay(self,'snowpackadvanced_variant','variant selection (includes a choice of specific models, DEFAULT, ANTARCTICA and JAPAN )'); % use 320 kg m-3 for fixed density
			fielddisplay(self,'snowpackadvanced_hn_density',{'Fixed value to be used as new snow density if a constant density model is chosen, otherwise the choices are "PARAMETERIZED" "EVENT" "MEASURED"'});
			% }}}
			disp(sprintf('\n	%s','General parameters:')); % {{{
			fielddisplay(self,'general_pluginpath','');
			fielddisplay(self,'general_buff_chunk_size','Size in days of a chunk of data to read at once.');
			fielddisplay(self,'general_buff_before','Alternate way of buffer centering: When rebuffering, the new date will be located BUFF_BEFORE days from the beginning of the buffer (therefore, it takes a value in days). ');
			% }}}
			disp(sprintf('\n	%s','Input  parameter:')); % {{{
			fielddisplay(self,'input_coordsys','coordinates in the Swiss Grid (http://geomatics.ladetto.ch/ch1903_wgs84_de.pdf). One of CH1903,UTM,UPS,PROJ4 or LOCAL');
			fielddisplay(self,'input_coordparam',' ');
			fielddisplay(self,'input_time_zone',' ');
			fielddisplay(self,'input_meteo','plugin for METEO data (one of BORMA,COSMO,GEOTOP,GRIB,GS,IMIS,SMET,SNOWPACK');
			fielddisplay(self,'input_meteopath','string containing the path to the xml files.');
			fielddisplay(self,'input_station1','Meteorology file for station number #');
			fielddisplay(self,'input_snowfile1','File name for the initial snow profile for station number #');
			% }}}
			disp(sprintf('\n	%s','Output parameters:')); % {{{
			fielddisplay(self,'output_coordsys','Coordinates in the Swiss Grid http://geomatics.ladetto.ch/ch1903_wgs84_de.pdf. One of CH1903,UTM,UPS,PROJ4 or LOCAL ');
			fielddisplay(self,'output_coordparam','');
			fielddisplay(self,'output_time_zone','');
			fielddisplay(self,'output_meteopath','Path to the outputs (this path MUST exist, it won''t be created)');
			fielddisplay(self,'output_experiment','Option to give an additional simulation specific output name to the run in addition to "STATION_NAME"');
			fielddisplay(self,'output_ts_write','Write meteo data out? (0 or 1)');
			fielddisplay(self,'output_ts_start','When to start writing meteo data out (offset, in days)');
			fielddisplay(self,'output_ts_days_between','How often to write meteo data out (in days: 3 hours=.125, 1 hour=4.1666e-2)');
			fielddisplay(self,'output_profile','How to write the profiles (default: ASCII, choice is ASCII,IMIS or ASCII IMIS)');
)');
			fielddisplay(self,'output_prof_write','Write profile data out? (0 or 1) ');
			fielddisplay(self,'output_prof_start','When to start writing profile data out (offset, in days)');
			fielddisplay(self,'output_prof_days_between','How often to write profile data out (in days: 3 hours=.125, 1 hour=4.1666e-2)');
			% }}}
			disp(sprintf('\n	%s','Interpolations1d parameters:')); % {{{
			fielddisplay(self,'interpolations1d_window_size','Affects resampling: expresses (in seconds) how far a valid point can be searched for when re-interpolating a missing value'); 
			fielddisplay(self,'interpolations1d_hnw_resample','NONE, NEAREST_NEIGHBOUR, ACCUMULATE or LINEAR');
 ');
			fielddisplay(self,'interpolations1d_hs_resample','Mean average processing. The mean average filter returns the mean value of all values within a user given time window. (NONE, NEAREST_NEIGHBOUR, ACCUMULATE or LINEAR)');
			fielddisplay(self,'interpolations1d_tsg_resample','Mean average processing. The mean average filter returns the mean value of all values within a user given time window.(NONE, NEAREST_NEIGHBOUR, ACCUMULATE or LINEAR)');
			fielddisplay(self,'interpolations1d_rho_hn_resample','(NONE, NEAREST_NEIGHBOUR, ACCUMULATE or LINEAR)');
			fielddisplay(self,'interpolations1d_vw_resample','(NONE, NEAREST_NEIGHBOUR, ACCUMULATE or LINEAR)');
			fielddisplay(self,'interpolations1d_vw_args','default nothing, otherwise, ''extrapolcate''');
			% }}}
			disp(sprintf('\n	%s','Filters parameters:')); % {{{
			fielddisplay(self,'filters_ta_filter1',' ');
			fielddisplay(self,'filters_ta_arg1','');
			fielddisplay(self,'filters_rh_filter1',' ');
			fielddisplay(self,'filters_rh_arg1','');
			fielddisplay(self,'filters_rh_filter2',' ');
			fielddisplay(self,'filters_rh_arg2','');
			fielddisplay(self,'filters_iswr_filter1',' ');
			fielddisplay(self,'filters_iswr_arg1','');
			fielddisplay(self,'filters_iswr_filter2',' ');
			fielddisplay(self,'filters_iswr_arg2','');
			fielddisplay(self,'filters_rswr_filter1',' ');
			fielddisplay(self,'filters_rswr_arg1','');
			fielddisplay(self,'filters_rswr_filter2',' ');
			fielddisplay(self,'filters_rswr_arg2','');

			%for ta between 190 and 280 k;
			fielddisplay(self,'filters_ilwr_filter1',' ');
			fielddisplay(self,'filters_ilwr_arg1','');
			fielddisplay(self,'filters_ilwr_filter2',' ');
			fielddisplay(self,'filters_ilwr_arg2','');
			fielddisplay(self,'filters_tss_filter1',' ');
			fielddisplay(self,'filters_tss_arg1','');
			fielddisplay(self,'filters_tsg_filter1',' ');
			fielddisplay(self,'filters_tsg_arg1','');
			fielddisplay(self,'filters_vw_filter1',' ');
			fielddisplay(self,'filters_vw_arg1','');
			fielddisplay(self,'filters_vw_filter2',' ');
			fielddisplay(self,'filters_vw_arg2','');
			% }}}

		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'object',self,'class','snowpack','fieldname','spcvx','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','snowpack','fieldname','spcvy','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','snowpack','fieldname','spcvz','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','snowpack','fieldname','restol','format','Double');
			WriteData(fid,prefix,'object',self,'class','snowpack','fieldname','reltol','format','Double');
			WriteData(fid,prefix,'object',self,'class','snowpack','fieldname','abstol','format','Double');
			WriteData(fid,prefix,'object',self,'class','snowpack','fieldname','isnewton','format','Integer');
			WriteData(fid,prefix,'object',self,'class','snowpack','fieldname','FSreconditioning','format','Double');
			WriteData(fid,prefix,'object',self,'class','snowpack','fieldname','maxiter','format','Integer');
			WriteData(fid,prefix,'object',self,'class','snowpack','fieldname','shelf_dampening','format','Integer');
			WriteData(fid,prefix,'object',self,'class','snowpack','fieldname','vertex_pairing','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',self,'class','snowpack','fieldname','penalty_factor','format','Double');
			WriteData(fid,prefix,'object',self,'class','snowpack','fieldname','rift_penalty_lock','format','Integer');
			WriteData(fid,prefix,'object',self,'class','snowpack','fieldname','rift_penalty_threshold','format','Integer');
			WriteData(fid,prefix,'object',self,'class','snowpack','fieldname','referential','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','snowpack','fieldname','requested_outputs','format','StringArray');
		end % }}}
	end
end
