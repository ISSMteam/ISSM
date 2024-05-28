%QMU class definition
%
%   Usage:
%      qmu=qmu();

classdef qmu
	properties (SetAccess=public)
		isdakota                    = 0;
		output                      = 0;
		variables                   = struct();
		correlation_matrix          = [];
		responses                   = struct();
		method                      = struct();
		params                      = struct();
		statistics                  = qmustatistics();
		results                     = struct();
		numberofresponses           = 0;
		variabledescriptors         = {};
		variablepartitions          = {};
		variablepartitions_npart    = [];
		variablepartitions_nt    = [];
		responsedescriptors         = {};
		responsepartitions          = {};
		responsepartitions_npart    = [];
		mass_flux_profile_directory = NaN;
		mass_flux_profiles          = NaN;
		mass_flux_segments          = {};
		adjacency                   = NaN;
		vertex_weight               = NaN;
	end
	methods (Static)
		function self = loadobj(self) % {{{
			% This function is directly called by matlab when a model object is
			% loaded. If the input is a struct it is an old version of this class and
			% old fields must be recovered (make sure they are in the deprecated
			% model properties)

			if verLessThan('matlab','7.9'),
				disp('Warning: your matlab version is old and there is a risk that load does not work correctly');
				disp('         if the model is not loaded correctly, rename temporarily loadobj so that matlab does not use it');

				% This is a Matlab bug: all the fields of md have their default value
				% Example of error message:
				% Warning: Error loading an object of class 'model':
				% Undefined function or method 'exist' for input arguments of type 'cell'
				%
				% This has been fixed in MATLAB 7.9 (R2009b) and later versions
			end

			if isstruct(self)
				%disp('Recovering qmu from older version');
				objstruct = self;
				self = structtoobj(qmu(),objstruct);

				%2019 Dec 7th
				%if isfield(objstruct,'partition'),      self.vpartition     = objstruct.partition;       end;
			end

		end% }}}
	end
	methods
		function self = extrude(self,md) % {{{
		end % }}}
		function self = qmu(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~md.qmu.isdakota, return; end

			version=IssmConfig('_DAKOTA_VERSION_'); version=str2num(version(1:3));

			if version < 6,
				if md.qmu.params.evaluation_concurrency~=1,
					md = checkmessage(md,['concurrency should be set to 1 when running dakota in library mode']);
				end
			else
				if ~strcmpi(self.params.evaluation_scheduling,'master'),
					md = checkmessage(md,['evaluation_scheduling in qmu.params should be set to ''master''']);
				end
				if md.cluster.nprocs()<=1,
					md = checkmessage(md,['in parallel library mode, Dakota needs to run on at least 2 cpus, 1 cpu for the master, 1 cpu for the slave. Modify md.cluster.np accordingly.']);
				end

				if self.params.processors_per_evaluation<1,
					md = checkmessage(md,['in parallel library mode, Dakota needs to run at least one slave on one cpu (md.qmu.params.processors_per_evaluation >=1)!']);
				end
				if mod(md.cluster.nprocs()-1,self.params.processors_per_evaluation),
					%md = checkmessage(md,['in parallel library mode, the requirement is for md.cluster.np = md.qmu.params.processors_per_evaluation * number_of_slaves, where number_of_slaves will automatically be determined by Dakota. Modify md.cluster.np accordingly']);
				end
			end

			%go through variables and check for consistency: 
			fv=fieldnames(self.variables);
			for i=1:length(fv),
				self.variables.(fv{i}).checkconsistency(md,solution,analyses);
			end

			%go through variables, and check that we have normal uncertains first, then uniform uncertains 
			%and finally histogram_bin_uncertain. Indeed, dakota will order them this way, and when we send 
			%partitions for scaled variables, they better show up in the order dakota is feeding them to us 
			%in InputUpdateFromDakotax!
			fv=fieldnames(self.variables); classlist={};
			for i=1:length(fv),
				classlist{i}=class(self.variables.(fv{i}));
			end
			n=0; u=0; h=0;
			for i=1:length(classlist),
				if strcmpi(classlist{i},'normal_uncertain')
					if (u~=0 | h~=0),
						error('normal_uncertain variables should be declared before uniform and histogram_bin uncertain variables');
					else
						n=1;
					end
				end
				if strcmpi(classlist{i},'uniform_uncertain')
					if (h~=0),
						error('uniform_uncertain variables should be declared before histogram_bin uncertain variables');
					else
						u=1;
					end
				end
				if strcmpi(classlist{i},'histogram_bin_uncertain')
					h=1;
				end

			end


		end % }}}
		function disp(self) % {{{
			disp(sprintf('   qmu parameters:'));

			fielddisplay(self,'isdakota','is qmu analysis activated?');
			fielddisplay(self,'output','are we outputting ISSM results, default is 0');
			for i=1:numel(self.variables)
				disp(sprintf('         variables%s:  (arrays of each variable class)',...
					string_dim(self.variables,i)));
				fnames=fieldnames(self.variables(i));
				maxlen=0;
				for j=1:numel(fnames)
					maxlen=max(maxlen,length(fnames{j}));
				end

				for j=1:numel(fnames)
					disp(sprintf(['            %-' num2str(maxlen+1) 's:    [%ix%i]    ''%s'''],...
						fnames{j},size(self.variables.(fnames{j})),class(self.variables.(fnames{j}))));
				end
			end
			for i=1:numel(self.responses)
				disp(sprintf('         responses%s:  (arrays of each response class)',...
					string_dim(self.responses,i)));
				fnames=fieldnames(self.responses(i));
				maxlen=0;
				for j=1:numel(fnames)
					maxlen=max(maxlen,length(fnames{j}));
				end

				for j=1:numel(fnames)
					disp(sprintf(['            %-' num2str(maxlen+1) 's:    [%ix%i]    ''%s'''],...
						fnames{j},size(self.responses.(fnames{j})),class(self.responses.(fnames{j}))));
				end
			end
			fielddisplay(self,'numberofresponses','number of responses')
			for i=1:numel(self.method);
				if strcmp(class(self.method(i)),'dakota_method')
					disp(sprintf('            method%s :    ''%s''',...
						string_dim(self.method,i),self.method(i).method));
				end
			end
			for i=1:numel(self.params)
				disp(sprintf('         params%s:  (array of method-independent parameters)',...
					string_dim(self.params,i)));
				fnames=fieldnames(self.params(i));
				maxlen=0;
				for j=1:numel(fnames)
					maxlen=max(maxlen,length(fnames{j}));
				end

				for j=1:numel(fnames)
					disp(sprintf(['            %-' num2str(maxlen+1) 's: %s'],...
						fnames{j},any2str(self.params(i).(fnames{j}))));
				end
			end
			for i=1:numel(self.results)
				disp(sprintf('         results%s:  (information from dakota files)',...
					string_dim(self.results,i)));
				fnames=fieldnames(self.results(i));
				maxlen=0;
				for j=1:numel(fnames)
					maxlen=max(maxlen,length(fnames{j}));
				end

				for j=1:numel(fnames)
					disp(sprintf(['            %-' num2str(maxlen+1) 's:    [%ix%i]    ''%s'''],...
						fnames{j},size(self.results.(fnames{j})),class(self.results.(fnames{j}))));
				end
			end
			fielddisplay(self,'variablepartitions','');
			fielddisplay(self,'variablepartitions_npart','');
			fielddisplay(self,'variablepartitions_nt','');
			fielddisplay(self,'variabledescriptors','');
			fielddisplay(self,'responsedescriptors','');
			fielddisplay(self,'method','array of dakota_method class');
			fielddisplay(self,'mass_flux_profile_directory','directory for mass flux profiles');
			fielddisplay(self,'mass_flux_profiles','list of mass_flux profiles');
			fielddisplay(self,'mass_flux_segments','');
			fielddisplay(self,'adjacency','');
			fielddisplay(self,'vertex_weight','weight applied to each mesh vertex');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'object',self,'fieldname','isdakota','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','output','format','Boolean');
			if ~self.isdakota,
				WriteData(fid,prefix,'data',false,'name','md.qmu.mass_flux_segments_present','format','Boolean');
				return;
			end
			WriteData(fid,prefix,'data',self.method.params.samples,'name','md.qmu.method.params.samples','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','numberofresponses','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','variabledescriptors','format','StringArray');
			WriteData(fid,prefix,'object',self,'fieldname','variablepartitions','format','MatArray');
			WriteData(fid,prefix,'object',self,'fieldname','variablepartitions_npart','format','IntMat','mattype',3);
			WriteData(fid,prefix,'object',self,'fieldname','variablepartitions_nt','format','IntMat','mattype',3);
			WriteData(fid,prefix,'object',self,'fieldname','responsedescriptors','format','StringArray');
			WriteData(fid,prefix,'object',self,'fieldname','responsepartitions','format','MatArray');
			WriteData(fid,prefix,'object',self,'fieldname','responsepartitions_npart','format','IntMat','mattype',3);
			if ~isempty(self.mass_flux_segments),
				WriteData(fid,prefix,'data',self.mass_flux_segments,'name','md.qmu.mass_flux_segments','format','MatArray');
				flag=true;
			else
				flag=false;
			end
			WriteData(fid,prefix,'data',flag,'name','md.qmu.mass_flux_segments_present','format','Boolean');
			self.statistics.marshall(prefix,md,fid);

		end % }}}
		function savemodeljs(self,fid,modelname) % {{{

			if self.isdakota,
				error('qmu savemodeljs error message: not supported yet!');
			end

		end % }}}
	end
end
