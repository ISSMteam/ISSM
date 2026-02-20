%SMBpddSicopolis Class definition
%
%   Usage:
%      SMBpddSicopolis=SMBpddSicopolis();

classdef SMBpddSicopolis
	properties (SetAccess=public)

		precipitation         = NaN;
		monthlytemperatures   = NaN;
		temperature_anomaly   = NaN;
		precipitation_anomaly = NaN;
		smb_corr              = NaN;
		desfac                = 0;
		s0p                   = NaN;
		s0t                   = NaN;
		rlaps                 = 0;
		isfirnwarming         = 0;
		pdd_fac_ice           = 0;
		pdd_fac_snow          = 0;
		steps_per_step        = 1
		averaging             = 0
		requested_outputs     = {};
	end
	methods
		function self = SMBpddSicopolis(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = extrude(self,md) % {{{
			self.precipitation=project3d(md,'vector',self.precipitation,'type','node');
			self.monthlytemperatures=project3d(md,'vector',self.monthlytemperatures,'type','node');
			self.temperature_anomaly=project3d(md,'vector',self.temperature_anomaly,'type','node');
			self.precipitation_anomaly=project3d(md,'vector',self.precipitation_anomaly,'type','node');
			self.smb_corr=project3d(md,'vector',self.smb_corr,'type','node');
			self.s0p=project3d(md,'vector',self.s0p,'type','node');
			self.s0t=project3d(md,'vector',self.s0t,'type','node');

		end % }}}
		function list = defaultoutputs(self,md) % {{{
			list = {'SmbMassBalance'};
		end % }}}
		function self = initialize(self,md) % {{{

			if isnan(self.s0p),
				self.s0p=zeros(md.mesh.numberofvertices,1);
				disp('      no SMBpddSicopolis.s0p specified: values set as zero');
			end
			if isnan(self.s0t),
				self.s0t=zeros(md.mesh.numberofvertices,1);
				disp('      no SMBpddSicopolis.s0t specified: values set as zero');
			end
			if isnan(self.temperature_anomaly),
				self.temperature_anomaly=zeros(md.mesh.numberofvertices,1);
				disp('      no SMBpddSicopolis.temperature_anomaly specified: values set as zero');
			end
			if isnan(self.precipitation_anomaly),
				self.precipitation_anomaly=ones(md.mesh.numberofvertices,1);
				disp('      no SMBpddSicopolis.precipitation_anomaly specified: values set as ones');
			end
			if isnan(self.smb_corr),
				self.smb_corr=zeros(md.mesh.numberofvertices,1);
				disp('      no SMBpddSicopolis.smb_corr specified: values set as zero');
			end

		end % }}}
		function self = setdefaultparameters(self) % {{{

			self.isfirnwarming = 1;
			self.desfac        = -log(2.0)/1000;
			self.rlaps         = 7.4;
			self.pdd_fac_ice   = 7.28;
			self.pdd_fac_snow  = 2.73;
			self.requested_outputs={'default'};

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if (strcmp(solution,'TransientSolution') & md.transient.issmb == 0), return; end

			if ismember('MasstransportAnalysis',analyses),
				md = checkfield(md,'fieldname','smb.desfac','<=',1,'numel',1);
				md = checkfield(md,'fieldname','smb.s0p','>=',0,'NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
				md = checkfield(md,'fieldname','smb.s0t','>=',0,'NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
				md = checkfield(md,'fieldname','smb.rlaps','>=',0,'numel',1);
				md = checkfield(md,'fieldname','smb.monthlytemperatures','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 12]);
				md = checkfield(md,'fieldname','smb.precipitation','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 12]);
				md = checkfield(md,'fieldname','smb.pdd_fac_ice','>',0,'numel',1);
				md = checkfield(md,'fieldname','smb.pdd_fac_snow','>',0,'numel',1);

			end
			md = checkfield(md,'fieldname','smb.steps_per_step','>=',1,'numel',[1]);
			md = checkfield(md,'fieldname','smb.averaging', 'numel', [1], 'values', [0, 1, 2]);
			md = checkfield(md,'fieldname','smb.requested_outputs','stringrow',1);

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   surface forcings parameters:'));

			disp(sprintf('\n   SICOPOLIS PDD scheme (Calov & Greve, 2005) :'));
			fielddisplay(self,'monthlytemperatures','monthly surface temperatures [K]');
			fielddisplay(self,'precipitation','monthly surface precipitation [m/yr water eq]');
			fielddisplay(self,'temperature_anomaly','anomaly to monthly reference temperature (additive; [K])');
			fielddisplay(self,'precipitation_anomaly','anomaly to monthly precipitation (multiplicative, e.g. q=q0*exp(0.070458*DeltaT) after Huybrechts (2002)); [no unit])');
			fielddisplay(self,'smb_corr','correction of smb after PDD call [m/a]');
			fielddisplay(self,'s0p','should be set to elevation from precip source (between 0 and a few 1000s m, default is 0) [m]');
			fielddisplay(self,'s0t','should be set to elevation from temperature source (between 0 and a few 1000s m, default is 0) [m]');
			fielddisplay(self,'rlaps','present day lapse rate (default is 7.4 degree/km)');
			fielddisplay(self,'desfac','desertification elevation factor (default is -log(2.0)/1000)');
			fielddisplay(self,'isfirnwarming','is firnwarming (Reeh 1991) activated (0 or 1, default is 1)');
			fielddisplay(self, 'steps_per_step', 'number of smb steps per time step');
			fielddisplay(self,'averaging','averaging methods from short to long steps');
			fielddisplay(self,'pdd_fac_ice','Pdd factor for ice for all the domain [mm ice equiv/day/degree C]');
			fielddisplay(self,'pdd_fac_snow','Pdd factor for snow for all the domain [mm ice equiv/day/degree C]');
			disp(sprintf('%51s  0: Arithmetic (default)',' '));
			disp(sprintf('%51s  1: Geometric',' '));
			disp(sprintf('%51s  2: Harmonic',' '));
			fielddisplay(self,'requested_outputs','additional outputs requested (TemperaturePDD, SmbAccumulation, SmbMelt)');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.smb.model','data',10,'format','Integer');

			WriteData(fid,prefix,'object',self,'class','smb','fieldname','isfirnwarming','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','desfac','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','s0p','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','s0t','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','rlaps','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','pdd_fac_ice','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','pdd_fac_snow','format','Double');

			WriteData(fid,prefix,'object',self,'class','smb','fieldname','monthlytemperatures','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','precipitation','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','temperature_anomaly','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','precipitation_anomaly','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','smb_corr','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid, prefix, 'object', self, 'fieldname', 'steps_per_step', 'format', 'Integer');
			WriteData(fid, prefix, 'object', self, 'fieldname', 'averaging', 'format', 'Integer')

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
