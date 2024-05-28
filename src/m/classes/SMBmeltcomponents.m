%SMBmeltcomponents Class definition
%
%   Usage:
%      SMBmeltcomponents=SMBmeltcomponents();

classdef SMBmeltcomponents
	properties (SetAccess=public)
		accumulation = NaN;
		evaporation = NaN;
		melt = NaN;
		refreeze = NaN;
		steps_per_step = 1;
		averaging = 0;
		requested_outputs= {};
	end
	methods
		function self = SMBmeltcomponents(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   surface forcings parameters with melt (SMB=accumulation-evaporation-melt+refreeze) :'));
			fielddisplay(self,'accumulation','accumulated snow [m/yr ice eq]');
			fielddisplay(self,'evaporation','amount of ice lost to evaporative processes [m/yr ice eq]');
			fielddisplay(self,'melt','amount of ice melt in ice column [m/yr ice eq]');
			fielddisplay(self,'refreeze','amount of ice melt refrozen in ice column [m/yr ice eq]');
			fielddisplay(self, 'steps_per_step', 'number of smb steps per time step');
			fielddisplay(self, 'averaging', 'averaging methods from short to long steps');
			disp(sprintf('%51s  0: Arithmetic (default)',' '));
			disp(sprintf('%51s  1: Geometric',' '));
			disp(sprintf('%51s  2: Harmonic',' '));
			fielddisplay(self,'requested_outputs','additional outputs requested');
		end % }}}
		function self = extrude(self,md) % {{{

			self.accumulation=project3d(md,'vector',self.accumulation,'type','node');
			self.evaporation=project3d(md,'vector',self.evaporation,'type','node');
			self.melt=project3d(md,'vector',self.melt,'type','node');
			self.refreeze=project3d(md,'vector',self.refreeze,'type','node');

		end % }}}
		function list = defaultoutputs(self,md) % {{{
			list = {'SmbMassBalance'};
		end % }}}
		function self = initialize(self,md) % {{{

			if isnan(self.accumulation)
				self.accumulation=zeros(md.mesh.numberofvertices,1);
				disp('      no smb.accumulation specified: values set as zero');
			end
			if isnan(self.evaporation)
				self.evaporation=zeros(md.mesh.numberofvertices,1);
				disp('      no smb.evaporation specified: values set as zero');
			end
			if isnan(self.refreeze)
				self.refreeze=zeros(md.mesh.numberofvertices,1);
				disp('      no smb.refreeze specified: values set as zero');
			end
			if isnan(self.melt)
				self.melt=zeros(md.mesh.numberofvertices,1);
				disp('      no smb.melt specified: values set as zero');
			end

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ismember('MasstransportAnalysis',analyses),
				md = checkfield(md,'fieldname','smb.accumulation','timeseries',1,'NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','smb.evaporation','timeseries',1,'NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','smb.refreeze','timeseries',1,'NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','smb.melt','timeseries',1,'NaN',1,'Inf',1);
			end
			if ismember('BalancethicknessAnalysis',analyses),
				md = checkfield(md,'fieldname','smb.accumulation','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','smb.evaporation','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','smb.refreeze','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','smb.melt','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
			end
			md = checkfield(md,'fieldname','smb.steps_per_step','>=',1,'numel',[1]);
			md = checkfield(md,'fieldname','smb.averaging','numel',[1],'values',[0 1 2]);
			md = checkfield(md,'fieldname','smb.requested_outputs','stringrow',1);
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.smb.model','data',3,'format','Integer');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','accumulation','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','evaporation','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','melt','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','refreeze','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid, prefix, 'object', self, 'fieldname', 'steps_per_step', 'format', 'Integer');
			WriteData(fid, prefix, 'object', self, 'fieldname', 'averaging', 'format', 'Integer');

			%process requested outputs
			outputs = self.requested_outputs;
			pos  = find(ismember(outputs,'default'));
			if ~isempty(pos),
				outputs(pos) = [];                         %remove 'default' from outputs
				outputs      = [outputs defaultoutputs(self,md)]; %add defaults
			end
			WriteData(fid,prefix,'data',outputs,'name','md.smb.requested_outputs','format','StringArray');

		end % }}}
		function self = setdefaultparameters(self) % {{{

			%output default:
			self.requested_outputs={'default'};

		end % }}}
	end
end
