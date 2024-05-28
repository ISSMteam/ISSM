%SMBsemic Class definition
%
%   Usage:
%      SMBsemic=SMBsemic();

classdef SMBsemic
	properties (SetAccess=public)
		dailysnowfall		= NaN;
		dailyrainfall		= NaN;
		dailydsradiation	= NaN;
		dailydlradiation	= NaN;
		dailywindspeed		= NaN;
		dailypressure		= NaN;
		dailyairdensity   = NaN;
		dailyairhumidity	= NaN;
		dailytemperature	= NaN;

		Tamp              = NaN;
		mask              = NaN;
		hice              = NaN;
		hsnow             = NaN;
		qmr               = NaN;
		desfac				= 0;
		desfacElevation   = 0;
		rlaps				= 0;
		rdl					= 0;
		s0gcm				= NaN;
		steps_per_step		= 1;
		averaging			= 0;
		requested_outputs	= {};

		hcrit             = 0;
		rcrit             = 0;

		% albedo
		albedo            = 0; % required for first energy balance calculation of SEMIC
		albedo_snow       = 0; % required for ISBA method
		albedo_scheme     = 0; 
		alb_smax = NaN;
		alb_smin = NaN;
		albi = NaN;
		albl = NaN;

		% albedo parameters depending on albedo_scheme
		% for slater 
		tmin = NaN;
		tmax = NaN;

		% for isba & denby method
		mcrit = NaN;

		% for isba
		tau_a = NaN;
		tau_f = NaN;
		wcrit = NaN;

		% for alex
		tmid = NaN;
		afac = NaN;

		% method
		ismethod  = 0;
	end
	methods
		function self = SMBsemic(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = extrude(self,md) % {{{
			self.dailysnowfall=project3d(md,'vector',self.dailysnowfall,'type','node');
			self.dailyrainfall=project3d(md,'vector',self.dailyrainfall,'type','node');
			self.dailydsradiation=project3d(md,'vector',self.dailydsradiation,'type','node');
			self.dailydlradiation=project3d(md,'vector',self.dailydlradiation,'type','node');
			self.dailywindspeed=project3d(md,'vector',self.dailywindspeed,'type','node');
			self.dailypressure=project3d(md,'vector',self.dailypressure,'type','node');
			self.dailyairdensity=project3d(md,'vector',self.dailyairdensity,'type','node');
			self.dailyairhumidity=project3d(md,'vector',self.dailyairhumidity,'type','node');
			self.dailytemperature=project3d(md,'vector',self.dailytemperature,'type','node');
			self.s0gcm=project3d(md,'vector',self.s0gcm,'type','node');

		end % }}}
		function list = defaultoutputs(self,md) % {{{
			list = {'SmbMassBalance'};
		end % }}}
		function list = outputlists(self,md) % {{{
			if self.ismethod == 1
				list = {'default','SmbMassBalance','SmbMassBalanceSnow','SmbMassBalanceIce',...
					'SmbMelt','SmbRefreeze','SmbAccumulation',...
					'SmbHIce','SmbHSnow','SmbAlbedo','SmbAlbedoSnow','TemperatureSEMIC',...
					'SmbSemicQmr','TotalSmb','TotalSmbMelt','TotalSmbRefreeze'};
			else
				list = {'default','SmbMassBalance'};
			end
		end % }}}
		function self = initialize(self,md) % {{{
			% Explain
			%  initialize SEMIC smb values, such as s0gcm(surface elevation), albedo,
			% albedo_snow, hice, hsnow, Tamp... values.
			% 
			%
			% Usage
			%  md.smb = initialize(md.smb,md);

			if isnan(self.s0gcm)
				if ~isnan(md.geometry.surface) & (numel(md.geometry.surface) == md.mesh.numberofvertices)
					self.s0gcm=md.geometry.surface;
					disp('      no SMBsemic.s0gcm specified: values from md.geometry.surface');
				else
					self.s0gcm=zeros(md.mesh.numberofvertices,1);
					disp('      no SMBsemic.s0gcm specified: values set as zero');
				end
			end
			if isnan(self.mask)
				self.mask = 2*ones(md.mesh.numberofvertices,1);
				disp('      no SMBsemic.mask specified: values set as 2 for ice');
			end

			% update each values.
			if isnan(self.Tamp)
				self.Tamp= 3*ones(md.mesh.numberofvertices,1);
				disp('      no SMBsemic.Tamp specified: values set as 3.0');
			end
			self.albedo     = 0.8*ones(md.mesh.numberofvertices,1);
			self.albedo_snow= 0.5*ones(md.mesh.numberofvertices,1);
			self.hice       = 10*ones(md.mesh.numberofvertices,1);
			self.hsnow      = 5*ones(md.mesh.numberofvertices,1);
			self.qmr        = zeros(md.mesh.numberofvertices,1);
		end % }}}
		function self = setdefaultparameters(self) % {{{

			% albedo parameters
			self.albedo_scheme   = 0;
			self.alb_smax = 0.79;
			self.alb_smin = 0.6;
			self.albi = 0.41;
			self.albl = 0.07;

			% albedo parameters for?
			% for slater
			self.tmin  = 263.15;
			self.tmax  = 273.15;
			% for isba & denby
			self.mcrit = 6e-8;
			% for isba
			self.tau_a = 0.008;
			self.tau_f = 0.24;
			self.wcrit = 15.0;
			% for alex
			self.tmid  = 273.35;
			self.afac  = -0.18;

			self.hcrit = 0.028;% from Krapp et al. (2017)
			self.rcrit = 0.85; % from Krapp et al. (2017)
		
			self.desfac		      = -log(2.0)/1000;
			self.desfacElevation = 2000;
			self.rlaps		      = 7.4;
			self.rdl			      = 0.29;

			self.ismethod        = 0;
			self.requested_outputs={'default'};
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ismember('MasstransportAnalysis',analyses)
				md = checkfield(md,'fieldname','smb.desfac','<=',1,'numel',1);
				md = checkfield(md,'fieldname','smb.s0gcm','>=',0,'NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
				md = checkfield(md,'fieldname','smb.rlaps','>=',0,'numel',1);
				md = checkfield(md,'fieldname','smb.rdl','>=',0,'numel',1);
				md = checkfield(md,'fieldname','smb.dailysnowfall','timeseries',1,'NaN',1,'Inf',1,'>=',0);
				md = checkfield(md,'fieldname','smb.dailyrainfall','timeseries',1,'NaN',1,'Inf',1,'>=',0);
				md = checkfield(md,'fieldname','smb.dailydsradiation','timeseries',1,'NaN',1,'Inf',1,'>=',0);
				md = checkfield(md,'fieldname','smb.dailydlradiation','timeseries',1,'NaN',1,'Inf',1,'>=',0);
				md = checkfield(md,'fieldname','smb.dailywindspeed','timeseries',1,'NaN',1,'Inf',1,'>=',0);
				md = checkfield(md,'fieldname','smb.dailypressure','timeseries',1,'NaN',1,'Inf',1,'>=',0);
				md = checkfield(md,'fieldname','smb.dailyairdensity','timeseries',1,'NaN',1,'Inf',1,'>=',0);
				md = checkfield(md,'fieldname','smb.dailyairhumidity','timeseries',1,'NaN',1,'Inf',1,'>=',0);
				md = checkfield(md,'fieldname','smb.dailytemperature','timeseries',1,'NaN',1,'Inf',1,'>=',0);

				% TODO: transient model should be merged with SEMIC developed by Ruckamp et al. (2018)

				md = checkfield(md,'fieldname','smb.ismethod','numel',1,'values',[0,1]);
				if self.ismethod == 1 % transient mode
					md = checkfield(md,'fieldname','smb.desfacElevation','>=',0,'numel',1);

					md = checkfield(md,'fieldname','smb.albedo_scheme','NaN',1,'Inf',1,'numel',1,'values',[0,1,2,3,4]);
					md = checkfield(md,'fieldname','smb.alb_smax','>=',0,'NaN',1,'Inf',1,'numel',1);
					md = checkfield(md,'fieldname','smb.mask','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1],'values',[0,1,2]);

					% initial values
					md = checkfield(md,'fieldname','smb.albedo','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1]);
					md = checkfield(md,'fieldname','smb.albedo_snow','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1]);
					md = checkfield(md,'fieldname','smb.alb_smax','>=',0,'<=',1,'NaN',1,'Inf',1,'numel',1);
					md = checkfield(md,'fieldname','smb.alb_smin','<=',1,'NaN',1,'Inf',1,'numel',1);
					md = checkfield(md,'fieldname','smb.albi','>=',0,'<=',1,'NaN',1,'Inf',1,'numel',1);
					md = checkfield(md,'fieldname','smb.albl','>=',0,'<=',1,'NaN',1,'Inf',1,'numel',1);
					md = checkfield(md,'fieldname','smb.hice','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1]);
					md = checkfield(md,'fieldname','smb.hsnow','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1]);
					md = checkfield(md,'fieldname','smb.qmr','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1]);
				end
			end
			md = checkfield(md,'fieldname','smb.steps_per_step','>=',1,'numel',[1]);
			md = checkfield(md,'fieldname','smb.averaging','numel',[1],'values',[0 1 2]);
			md = checkfield(md,'fieldname','smb.requested_outputs','stringrow',1);
			% check requested_outputs
			if self.ismethod==1
				for i = 1:length(self.requested_outputs)
					if ~any(strcmpi(self.requested_outputs{i},self.outputlists))
						error(sprintf('ERROR: %s requested_output is not available',self.requested_outputs{i}));
					end
				end
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   surface forcings parameters:'));

			disp(sprintf('   Interface for coupling GCM data to the energy balance model SEMIC (Krapp et al (2017) https://doi.org/10.5194/tc-11-1519-2017).'));
			disp(sprintf('   The implemented coupling uses daily mean GCM input to calculate yearly mean smb, accumulation, ablation, and surface temperature.'));
			disp(sprintf('   smb and temperatures are updated every year'));
			disp(sprintf('\n   SEMIC parameters:'));
			fielddisplay(self,'dailysnowfall','daily surface dailysnowfall [m/s]');
			fielddisplay(self,'dailyrainfall','daily surface dailyrainfall [m/s]');
			fielddisplay(self,'dailydsradiation','daily downwelling shortwave radiation [W/m2]');
			fielddisplay(self,'dailydlradiation','daily downwelling longwave radiation [W/m2]');
			fielddisplay(self,'dailywindspeed','daily surface wind speed [m/s]');
			fielddisplay(self,'dailypressure','daily surface pressure [Pa]');
			fielddisplay(self,'dailyairdensity','daily air density [kg/m3]');
			fielddisplay(self,'dailyairhumidity','daily air specific humidity [kg/kg]');
			fielddisplay(self,'dailytemperature','daily surface air temperature [K]');
			fielddisplay(self,'rlaps','present day lapse rate (default is 7.4 [degree/km]; Erokhina et al. 2017)');
			fielddisplay(self,'desfac','desertification elevation factor (default is -log(2.0)/1000 [1/km]; Vizcaino et al. 2010)');
			fielddisplay(self,'rdl','longwave downward radiation decrease (default is 0.29 [W/m^2/km]; Marty et al. 2002)');
			fielddisplay(self,'s0gcm','GCM reference elevation; (default is 0) [m]');
			fielddisplay(self,'albedo_scheme','albedom scheme. 0: none, 1: (default is 0)');

			fielddisplay(self,'ismethod','method for calculating SMB with SEMIC. Default version of SEMIC is really slow. 0: steady, 1: transient (default: 0)');
			if self.ismethod == 1 % transient mode
				fielddisplay(self,'desfacElevation','desertification elevation (default is 2000 m; Vizcaino et al. 2010)');
				fielddisplay(self,'Tamp','amplitude of diurnal cycle [K]');
				fielddisplay(self,'albedo','initial albedo [no unit]');
				fielddisplay(self,'albedo_snow','initial albedo for snow [no unit]');
				fielddisplay(self,'hice','initial thickness of ice [unit: m]');
				fielddisplay(self,'hsnow','initial thickness of snow [unit: m]');
				fielddisplay(self,'mask','masking for albedo. 0: ocean, 1: land, 2: ice (default: 2)');
				fielddisplay(self,'qmr','initial net energy difference between melt and refreeze in SEMIC [unit: W m^{-2}]. This variable can be set with zeros because net energy difference between melt and refreeze is dissipated fast.');
				fielddisplay(self,'hcrit','critical snow height for albedo [unit: m]');
				fielddisplay(self,'rcrit','critical refreezing height for albedo [no unit]');

				disp(sprintf('\nSEMIC albedo parameters.'));
				fielddisplay(self,'albedo_scheme','albedo scheme for SEMIC. 0: none, 1: slater, 2: isba, 3: denby, 4: alex (default is 0)');
				fielddisplay(self,'alb_smax','maximum snow albedo (default: 0.79)');
				fielddisplay(self,'alb_smin','minimum snow albedo (default: 0.6)');
				fielddisplay(self,'albi','background albedo for bare ice (default: 0.41)');
				fielddisplay(self,'albl','background albedo for bare land (default: 0.07)');
			end
			% albedo_scheme - 0: none, 1: slater, 2: isba, 3: denby, 4: alex.
         if self.albedo_scheme == 0
            disp(sprintf('\n\tSEMIC snow albedo parameter of None.'));
				disp(sprintf('\t   albedo of snow is updated from albedo snow max (alb_smax).'));
            disp(sprintf('\t   alb_snow = abl_smax'));
			elseif self.albedo_scheme == 1
				disp(sprintf('\n\tSEMIC snow albedo parameters of Slater et al, (1998).'));
				disp(sprintf('\t   alb = alb_smax - (alb_smax - alb_smin)*tm^(3.0)'))
				disp(sprintf('\t   tm  = 1 (tsurf > 273.15 K)'));
				disp(sprintf('\t         tm = f*(tsurf-tmin) (tmin <= tsurf < 273.15)'));
				disp(sprintf('\t         0 (tsurf < tmin)'));
				disp(sprintf('\t   f = 1/(273.15-tmin)'));
				fielddisplay(self,'tmin','minimum temperature for which albedo decline become effective. (default: 263.15 K)[unit: K])');
				fielddisplay(self,'tmax','maxmium temperature for which albedo decline become effective. This value should be fixed. (default: 273.15 K)[unit: K])');
			elseif self.albedo_scheme == 2
				disp(sprintf('\n\tSEMIC snow albedo parameters for ISBA.? where is citation?'));
				fielddisplay(self,'mcrit','critical melt rate (default: 6e-8) [unit: m/sec]');
				fielddisplay(self,'wcrit','critical liquid water content (default: 15) [unit: kg/m2]');
				fielddisplay(self,'tau_a','dry albedo decline [unit: 1/day]');
				fielddisplay(self,'tau_f','wet albedo decline [unit: 1/day]');
			elseif self.albedo_scheme == 3
				disp(sprintf('\n\tSEMIC snow albedo parameters for Denby et al. (2002 Tellus)'));
				fielddisplay(self,'mcrit','critical melt rate (default: 6e-8) [unit: m/sec]');
			elseif self.albedo_scheme == 4
				disp(sprintf('\n\tSEMIC snow albedo parameters for Alex.?'));
				fielddisplay(self,'afac','[unit: ?]');
				fielddisplay(self,'tmid','[unit: ?]');
			else
				error(sprintf('ERROR: %d is not supported albedom scheme.',self.albedo_scheme))
			end

			fielddisplay(self, 'steps_per_step', 'number of smb steps per time step');
			fielddisplay(self, 'averaging', 'averaging methods from short to long steps');
			disp(sprintf('%51s  0: Arithmetic (default)',' '));
			disp(sprintf('%51s  1: Geometric',' '));
			disp(sprintf('%51s  2: Harmonic',' '));
			fielddisplay(self,'requested_outputs','additional outputs requested');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.smb.model','data',12,'format','Integer');

			WriteData(fid,prefix,'object',self,'class','smb','fieldname','ismethod','format','Integer','values',[0, 1]);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','desfac','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','desfacElevation','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','s0gcm','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','rlaps','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','rdl','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','dailysnowfall','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','dailyrainfall','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','dailydsradiation','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','dailydlradiation','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','dailywindspeed','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','dailypressure','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','dailyairdensity','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','dailyairhumidity','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','dailytemperature','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
			% TODO: transient mode should be merged with SEMIC developed by Ruckamp et al. (2018).
			if self.ismethod == 1% transient mode
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','Tamp','format','DoubleMat','mattype',1);
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','mask','format','DoubleMat','mattype',1);
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','hice','format','DoubleMat','mattype',1);
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','hsnow','format','DoubleMat','mattype',1);
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','qmr','format','DoubleMat','mattype',1);

				WriteData(fid,prefix,'object',self,'class','smb','fieldname','hcrit','format','Double');
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','rcrit','format','Double');

				%albedo
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','albedo','format','DoubleMat','mattype',1);
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','albedo_snow','format','DoubleMat','mattype',1);
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','albedo_scheme','format','Integer');
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','albi','format','Double');
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','albl','format','Double');
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','alb_smin','format','Double');
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','alb_smax','format','Double');

				%albedo parameters for ?
				%for slater
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','tmin','format','Double');
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','tmax','format','Double');
				%for isba & denby
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','mcrit','format','Double');
				%for isba
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','wcrit','format','Double');
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','tau_a','format','Double');
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','tau_f','format','Double');
				%for alex
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','tmid','format','Double');
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','afac','format','Double');
			end

			WriteData(fid,prefix,'object',self,'fieldname','steps_per_step','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','averaging','format','Integer');

			%process requested outputs
			outputs = self.requested_outputs;
			pos  = find(ismember(outputs,'default'));
			if ~isempty(pos)
				outputs(pos) = []; %remove 'default' from outputs
				outputs      = [outputs defaultoutputs(self,md)]; %add defaults
			end
			WriteData(fid,prefix,'data',outputs,'name','md.smb.requested_outputs','format','StringArray');

		end % }}}
	end
end
