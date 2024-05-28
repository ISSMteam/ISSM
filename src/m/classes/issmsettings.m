%ISSMSETTINGS class definition
%
%   Usage:
%      issmsettings=issmsettings();

classdef issmsettings
	properties (SetAccess=public) 
		results_on_nodes         = {};
		io_gather                = 0;
		lowmem                   = 0;
		output_frequency         = 0;
		sb_coupling_frequency    = 0;
		checkpoint_frequency     = 0;
		waitonlock               = 0;
		upload_server            = '';
		upload_path              = '';
		upload_login             = '';
		upload_port              = 0;
		upload_filename          = '';
		solver_residue_threshold = 0;
	end
	methods (Static)
		function self = loadobj(self) % {{{
			% This function is directly called by matlab when a model object is
			% loaded. Update old properties here

			%2020 Oct 6
			if isstruct(self)
				objstruct = self;
				self = structtoobj(issmsettings(),objstruct);
				if isfield(objstruct,'recording_frequency')
					self.checkpoint_frequency = objstruct.recording_frequency;
				end
			end
		end % }}}
	end
	methods
		function self = issmsettings(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

			%are we short in memory ? (0 faster but requires more memory)
			self.lowmem=0;

			%i/o:
			self.io_gather=1;

			%results frequency by default every step
			self.output_frequency=1;

			%coupling frequency of the stress balance solver by default every step
			self.sb_coupling_frequency=1;
			
			%checkpoints frequency, by default never: 
			self.checkpoint_frequency=0;

			%this option can be activated to load automatically the results
			%onto the model after a parallel run by waiting for the lock file
			%N minutes that is generated once the solution has converged
			%0 to deactivate
			self.waitonlock=Inf;

			%upload options: 
			self.upload_port         = 0;

			%throw an error if solver residue exceeds this value
			self.solver_residue_threshold = 1e-6;

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			md = checkfield(md,'fieldname','settings.results_on_nodes','stringrow',1);
			md = checkfield(md,'fieldname','settings.io_gather','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','settings.lowmem','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','settings.output_frequency','numel',[1],'>=',1);
			md = checkfield(md,'fieldname','settings.sb_coupling_frequency','numel',[1],'>=',1);
			md = checkfield(md,'fieldname','settings.checkpoint_frequency','numel',[1],'>=',0);
			md = checkfield(md,'fieldname','settings.waitonlock','numel',[1]);
			md = checkfield(md,'fieldname','settings.solver_residue_threshold','numel',[1],'>',0);

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   general issmsettings parameters:'));

			fielddisplay(self,'results_on_nodes','list of output for which results will be output for all the nodes of each element, Use ''all'' for all output on nodes.');
			fielddisplay(self,'io_gather','I/O gathering strategy for result outputs (default 1)');
			fielddisplay(self,'lowmem','is the memory limited ? (0 or 1)');
			fielddisplay(self,'output_frequency','number of time steps between two saves (e.g., 5 means that results are only saved every 5 time steps)');
			fielddisplay(self,'sb_coupling_frequency','frequency at which StressBalance solver is coupled (default 1)');
			fielddisplay(self,'checkpoint_frequency','frequency at which the runs are being recorded, allowing for a restart');
			fielddisplay(self,'waitonlock','maximum number of minutes to wait for batch results (NaN to deactivate)');
			fielddisplay(self,'upload_server','server hostname where model should be uploaded');
			fielddisplay(self,'upload_path','path on server where model should be uploaded');
			fielddisplay(self,'upload_login','server login');
			fielddisplay(self,'upload_port','port login (default is 0)');
			fielddisplay(self,'upload_filename','unique id generated when uploading the file to server');
			fielddisplay(self,'solver_residue_threshold','throw an error if solver residue exceeds this value (NaN to deactivate)');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'data',self.results_on_nodes,'name','md.settings.results_on_nodes','format','StringArray');
			WriteData(fid,prefix,'object',self,'class','settings','fieldname','io_gather','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','settings','fieldname','lowmem','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','settings','fieldname','output_frequency','format','Integer');
			WriteData(fid,prefix,'object',self,'class','settings','fieldname','sb_coupling_frequency','format','Integer');
			WriteData(fid,prefix,'object',self,'class','settings','fieldname','checkpoint_frequency','format','Integer');
			WriteData(fid,prefix,'object',self,'class','settings','fieldname','waitonlock','data',self.waitonlock>0,'format','Boolean');
			WriteData(fid,prefix,'object',self,'class','settings','fieldname','solver_residue_threshold','format','Double');
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejscellstring(fid,[modelname '.settings.results_on_nodes'],self.results_on_nodes);
			writejsdouble(fid,[modelname '.settings.io_gather'],self.io_gather);
			writejsdouble(fid,[modelname '.settings.lowmem'],self.lowmem);
			writejsdouble(fid,[modelname '.settings.output_frequency'],self.output_frequency);
			writejsdouble(fid,[modelname '.settings.sb_coupling_frequency'],self.sb_coupling_frequency);
			writejsdouble(fid,[modelname '.settings.checkpoint_frequency'],self.checkpoint_frequency);
			writejsdouble(fid,[modelname '.settings.waitonlock'],self.waitonlock);
			writejsstring(fid,[modelname '.settings.upload_server'],self.upload_server);
			writejsstring(fid,[modelname '.settings.upload_path'],self.upload_path);
			writejsstring(fid,[modelname '.settings.upload_login'],self.upload_login);
			writejsdouble(fid,[modelname '.settings.upload_port'],self.upload_port);
			writejsstring(fid,[modelname '.settings.upload_filename'],self.upload_filename);
			writejsstring(fid,[modelname '.settings.solver_residue_threshold'],self.solver_residue_threshold);
		end % }}}
	end
end
