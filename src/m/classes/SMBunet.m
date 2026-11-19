%SMBunet class definition
%
%   Usage:
%      SMBunet=SMBunet();

classdef SMBunet
	properties (SetAccess=public)

		x                     = NaN;
		y                     = NaN;
		steps_per_step        = 1;
		requested_outputs     = {};
		averaging             = 0;
      module_dir = '';
		pt_name = '';
		py_name = '';
	end
	methods
		function self = SMBunet(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = extrude(self,md) % {{{
			disp('     projection is not implemented.');
		end % }}}
		function list = defaultoutputs(self,md) % {{{
			list = {};
		end % }}}
		function self = initialize(self,md) % {{{
			if isnan(self.x)
				self.x=zeros(100,1);
				disp('      no smb.x specified: values set as zero');
			end
			if isnan(self.y)
				self.y=zeros(100,1);
				disp('      no smb.y specified: values set as zero');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{
			self.requested_outputs={'default'};
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			if (strcmp(solution,'TransientSolution') & md.transient.issmb == 0), return; end
			if ismember('MasstransportAnalysis',analyses)
			md = checkfield(md,'fieldname','smb.x','timeseries',0,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','smb.y','timeseries',0,'NaN',1,'Inf',1);
			end
			md = checkfield(md,'fieldname','smb.steps_per_step','>=',1,'numel',[1]);
			md = checkfield(md,'fieldname','smb.requested_outputs','stringrow',1);
			md = checkfield(md,'fieldname','smb.averaging','numel',[1],'values',[0 1 2]);
			md = checkfield(md,'fieldname','smb.module_dir','filepath',1);
			md = checkfield(md,'fieldname','smb.py_name','string',1);
			md = checkfield(md,'fieldname','smb.pt_name','string',1);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('\n  smb emulator based on U-Net :'));
         fielddisplay(self, 'x', 'U-Net grid x coordinates as a vector');
         fielddisplay(self, 'y', 'U-Net grid y coordinates as a vector');
			fielddisplay(self, 'steps_per_step', 'number of smb steps per time step');
			fielddisplay(self,'requested_outputs','additional outputs requested');
			fielddisplay(self,'averaging','averaging methods from short to long steps');
         fielddisplay(self,'module_dir', 'directory of the emulator module');
			fielddisplay(self,'pt_name', 'name of the checkpoint file for pre-trained ML model');
			fielddisplay(self,'py_name', 'name of the python file that defines ML architecture');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'name','md.smb.model','data',21,'format','Integer');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','x','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','y','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','module_dir','format','String')
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','pt_name','format','String');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','py_name','format','String');
			WriteData(fid,prefix,'object',self,'fieldname','steps_per_step','format','Integer');
			 WriteData(fid,prefix,'object',self,'fieldname','averaging','format','Integer');

			%process requested outputs
			outputs = self.requested_outputs;
			pos  = find(ismember(outputs,'default'));
			if ~isempty(pos)
				outputs(pos) = [];                         %remove 'default' from outputs
				outputs      = [outputs defaultoutputs(self,md)]; %add defaults
			end
			WriteData(fid,prefix,'data',outputs,'name','md.smb.requested_outputs','format','StringArray');

		end % }}}
	end
end
