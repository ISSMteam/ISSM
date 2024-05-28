%SMBgradientscomponents Class definition
%
%   Usage:
%      SMBgradientscomponents=SMBgradientscomponents();

classdef SMBgradientscomponents
	properties (SetAccess=public)

		accuref           = NaN;
		accualti          = NaN;
		accugrad          = NaN;
		runoffref         = NaN;
		runoffalti        = NaN;
		runoffgrad        = NaN;
		steps_per_step    = 1;
		averaging         = 0;
		requested_outputs = {};
	end
	methods
		function self=SMBgradientscomponents(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}

		function self = extrude(self,md) % {{{

		%Nothing for now

		end % }}}
		function list = defaultoutputs(self,md) % {{{
			list = {'SmbMassBalance'};
		end % }}}
		function self = initialize(self,md) % {{{

		%Nothing done for now

		end % }}}
		function self = setdefaultparameters(self) % {{{

			%output default:
			self.requested_outputs={'default'};

		end % }}}
		function	md=checkconsistency(self,md,solution,analyses) % {{{
			if ismember('MasstransportAnalysis',analyses),
				md = checkfield(md,'fieldname','smb.accuref','singletimeseries',1,'NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','smb.accualti','numel',1,'NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','smb.accugrad','singletimeseries',1,'NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','smb.runoffref','singletimeseries',1,'NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','smb.runoffalti','numel',1,'NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','smb.runoffgrad','singletimeseries',1,'NaN',1,'Inf',1);
			end
			md = checkfield(md,'fieldname','smb.steps_per_step','>=',1,'numel',[1]);
			md = checkfield(md,'fieldname','smb.averaging','numel',[1],'values',[0 1 2]);
			md = checkfield(md,'fieldname','masstransport.requested_outputs','stringrow',1);
		end % }}}

		function disp(self) % {{{
			disp(sprintf('   surface forcings parameters:'));
			disp(sprintf('\n   SMB gradients components parameters:'));
			fielddisplay(self,'accuref',' reference value of the accumulation');
			fielddisplay(self,'accualti',' Altitude at which the accumulation is equal to the reference value');
			fielddisplay(self,'accugrad',' Gradient of the variation of the accumulation (0 for uniform accumulation)');
			fielddisplay(self,'runoffref',' reference value of the runoff m w.e. y-1 (temperature times ddf)');
			fielddisplay(self,'runoffalti',' Altitude at which the runoff is equal to the reference value');
			fielddisplay(self,'runoffgrad',' Gradient of the variation of the runoff (0 for uniform runoff) m w.e. m-1 y-1 (lpase rate times ddf)');
			fielddisplay(self, 'steps_per_step', 'number of smb steps per time step');
			fielddisplay(self, 'averaging', 'averaging methods from short to long steps');
			disp(sprintf('%51s  0: Arithmetic (default)',' '));
			disp(sprintf('%51s  1: Geometric',' '));
			disp(sprintf('%51s  2: Harmonic',' '));
			fielddisplay(self,'requested_outputs','additional outputs requested');

		end % }}}
		function marshall(self,prefix,md,fid)    % {{{

			WriteData(fid,prefix,'name','md.smb.model','data',11,'format','Integer');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','accuref','format','DoubleMat','mattype',1,'timeserieslength',2,'yts',md.constants.yts,'scale',1./md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','accualti','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','accugrad','format','DoubleMat','mattype',1,'timeserieslength',2,'yts',md.constants.yts,'scale',1./md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','runoffref','format','DoubleMat','mattype',1,'timeserieslength',2,'yts',md.constants.yts,'scale',1./md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','runoffalti','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','runoffgrad','format','DoubleMat','mattype',1,'timeserieslength',2,'yts',md.constants.yts,'scale',1./md.constants.yts);
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
	end
end
