%SMBdebrisEvatt Class definition
%
%   Usage:
%      SMBdebrisEvatt=SMBdebrisEvatt();

classdef SMBdebrisEvatt
	properties (SetAccess=public)

		precipitation  = NaN;
		temperature   = NaN;
		dsradiation    = NaN;
                dlradiation    = NaN;
                windspeed      = NaN;
                airhumidity    = NaN;
		precipitation_anomaly = NaN;
		temperature_anomaly   = NaN;
		dsradiation_anomaly   = NaN;
		dlradiation_anomaly   = NaN;
		windspeed_anomaly     = NaN;
		airhumidity_anomaly   = NaN;
		s0t                   = NaN;
		snowheight            = NaN;
		qlaps                 = 0;
		rlaps                 = 0;
		dsgrad		      = 0;
		dlgrad		      = 0;
		windspeedgrad	      = 0;
		humiditygrad	      = 0;
		isAnderson	      = 0;
		iscryokarst	      = 0;
		AndersonD0	      = 0;
		steps_per_step        = 1;
		averaging             = 0;
		requested_outputs     = {};
		icealbedo	      = NaN;
		snowalbedo            = NaN;
		debrisalbedo	      = NaN;
	end
	methods
		function self = SMBdebrisEvatt(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = extrude(self,md) % {{{
			self.precipitation=project3d(md,'vector',self.precipitation,'type','node');
			self.temperature=project3d(md,'vector',self.temperature,'type','node');
			self.dsradiation=project3d(md,'vector',self.dsradiation,'type','node');
			self.dlradiation=project3d(md,'vector',self.dlradiation,'type','node');
			self.windspeed=project3d(md,'vector',self.windspeed,'type','node');
			self.airhumidity=project3d(md,'vector',self.airhumidity,'type','node');
			self.temperature_anomaly=project3d(md,'vector',self.temperature_anomaly,'type','node');
                        self.precipitation_anomaly=project3d(md,'vector',self.precipitation_anomaly,'type','node');
			self.dsradiation_anomaly=project3d(md,'vector',self.temperature_anomaly,'type','node');
	                self.dlradiation_anomaly=project3d(md,'vector',self.temperature_anomaly,'type','node');
         	        self.windspeed_anomaly=project3d(md,'vector',self.temperature_anomaly,'type','node');
                	self.airhumidity_anomaly=project3d(md,'vector',self.temperature_anomaly,'type','node');

			self.s0t=project3d(md,'vector',self.s0t,'type','node');
			self.snowheight=project3d(md,'vector',self.snowheight,'type','node');

		end % }}}
		function list = defaultoutputs(self,md) % {{{
			list = {'SmbMassBalance'};
		end % }}}
		function self = initialize(self,md) % {{{

			if isnan(self.s0t),
				self.s0t=zeros(md.mesh.numberofvertices,1);
				disp('      no SMBdebrisEvatt.s0t specified: values set as zero');
			end

		end % }}}
		function self = setdefaultparameters(self) % {{{

			self.qlaps         = 0.1;
			self.rlaps         = 7.4;
			self.dsgrad	   = 13.;
			self.dlgrad	   = 29;
			self.windspeedgrad = -0.2;
			self.humiditygrad  = 0;
			self.icealbedo	   = 0.3;
			self.snowalbedo    = 0.75;
		 	self.debrisalbedo  = 0.07;
			self.isAnderson    = 0;
			self.iscryokarst   = 0;
			self.AndersonD0    = 0.5;
			self.requested_outputs={'default'};

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if (strcmp(solution,'TransientSolution') & md.transient.issmb == 0), return; end

			if ismember('MasstransportAnalysis',analyses),
				md = checkfield(md,'fieldname','smb.temperature','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 365]);
				md = checkfield(md,'fieldname','smb.precipitation','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 365]);
				md = checkfield(md,'fieldname','smb.dsradiation','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 365]);
				md = checkfield(md,'fieldname','smb.dlradiation','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 365]);
				md = checkfield(md,'fieldname','smb.windspeed','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 365]);
				md = checkfield(md,'fieldname','smb.airhumidity','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 365]);
				md = checkfield(md,'fieldname','smb.snowheight','>=',0,'NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
			end
			md = checkfield(md,'fieldname','smb.steps_per_step','>=',1,'numel',[1]);
			md = checkfield(md,'fieldname','smb.averaging', 'numel', [1], 'values', [0, 1, 2]);
			md = checkfield(md,'fieldname','smb.requested_outputs','stringrow',1);
			md = checkfield(md,'fieldname','smb.icealbedo','>=',0,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','smb.snowalbedo','>=',0,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','smb.debrisalbedo','>=',0,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','smb.isAnderson','numel', [1], 'values', [0, 1]);
			md = checkfield(md,'fieldname','smb.iscryokarst','numel', [1], 'values', [0, 1]);

		end % }}}
		function disp(self) % {{{
	
			disp(sprintf('Evatt et al. (2015) Debris Model (doi: 10.3189/2015JoG14J235)'));
			disp(sprintf('If isAnderson==1 -> Eq. 6 from Ferguson & Vieli (2021) ist taken (https://doi.org/10.5194/tc-15-3377-2021)'));
			disp(sprintf('If iscryokarst==1 -> Eqs. 9,10 from Ferguson & Vieli (2021) are taken (https://doi.org/10.5194/tc-15-3377-2021)'));
			disp(sprintf('Clean-ice SMB is taken from the Evatt et al. (2015) EBM with debris=0'));

			fielddisplay(self,'isAnderson','do we use the Anderson parametrization (default is 0)');
			fielddisplay(self,'iscryokarst','do we use a cryokarst parametrization (default is 0)');
			fielddisplay(self,'temperature',' surface temperature [K]');
			fielddisplay(self,'precipitation',' surface precipitation [m/yr water eq]');
			fielddisplay(self,'dsradiation',' downwelling shortwave radiation [W m-2]');
                        fielddisplay(self,'dlradiation',' downwelling longwave radiation [W m-2]');
			fielddisplay(self,'windspeed',' surface wind speed [m s-1]');
                        fielddisplay(self,'airhumidity',' near-surface specific humidity [kg kg-1]');
			fielddisplay(self,'temperature_anomaly','anomaly to  reference temperature (additive)');
                        fielddisplay(self,'precipitation_anomaly','anomaly to  precipitation (multiplicative)');
                        fielddisplay(self,'dsradiation_anomaly','anomaly to  reference downwelling shortwave radiation');
                        fielddisplay(self,'dlradiation_anomaly','anomaly to  reference downwelling longwave radiation (additive');
                        fielddisplay(self,'windspeed_anomaly','anomaly to  reference surface wind speed (additive)');
                        fielddisplay(self,'airhumidity_anomaly','anomaly to  reference near-surface specific humidity (additive)');

			fielddisplay(self,'s0t','should be set to elevation from RCM/GCM source (between 0 and a few 1000s m, default is 0) [m]');
			fielddisplay(self,'snowheight','guess of snowheight at the end of the summer, will be further evolved');
			fielddisplay(self,'rlaps','present day temperature lapse rate (default is 7.4 degree/km)');
			fielddisplay(self,'dsgrad','present day SW height gradient (default is 1.3 W/m^2/km)');
			fielddisplay(self,'dlgrad','present day LW height gradient (default is 2.9 W/m^2/km)');
			fielddisplay(self,'windspeedgrad','present day wind speed height gradient (default is 0.02 m/s/km)');
			fielddisplay(self,'humiditygrad','present day humidity height gradient (default is 0)');
			fielddisplay(self,'qlaps','precip change (default is 0.1/km');
			fielddisplay(self,'icealbedo','albedo for ice (default is 0.3)');
			fielddisplay(self,'snowalbedo','albedo for snow (default is 0.75)');
			fielddisplay(self,'debrisalbedo','albedo for debris (default is 0.07)');
			fielddisplay(self,'AndersonD0','parameter to represent the debris effect (default is 0.5)');
			fielddisplay(self,'steps_per_step', 'number of smb steps per time step');
			fielddisplay(self,'averaging','averaging methods from short to long steps');
			disp(sprintf('%51s  0: Arithmetic (default)',' '));
			disp(sprintf('%51s  1: Geometric',' '));
			disp(sprintf('%51s  2: Harmonic',' '));
			fielddisplay(self,'requested_outputs','additional outputs requested (TemperaturePDD, SmbAccumulation, SmbMelt, SmbSummerMelt, SmbAlbedo, SmbSummerAlbedo, SmbSnowheight)');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.smb.model','data',14,'format','Integer');

			WriteData(fid,prefix,'object',self,'class','smb','fieldname','qlaps','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','s0t','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','snowheight','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','rlaps','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','dsgrad','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','dlgrad','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','windspeedgrad','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','humiditygrad','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','icealbedo','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','snowalbedo','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','debrisalbedo','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','isAnderson','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','iscryokarst','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','AndersonD0','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','temperature','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','precipitation','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','dsradiation','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','dlradiation','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','windspeed','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','airhumidity','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','temperature_anomaly','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
                        WriteData(fid,prefix,'object',self,'class','smb','fieldname','precipitation_anomaly','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','dsradiation_anomaly','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','dlradiation_anomaly','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','windspeed_anomaly','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','airhumidity_anomaly','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','steps_per_step','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','averaging','format','Integer')

			%process requested outputs
			outputs = self.requested_outputs;
			pos  = find(ismember(outputs,'default'));
			if ~isempty(pos),
				outputs(pos) = [];                         %remove 'default' from outputs
				outputs      = [outputs defaultoutputs(self,md)]; %add defaults
			end
			WriteData(fid,prefix,'data',outputs,'name','md.smb.requested_outputs','format','StringArray');

		end % }}}
	end
end
