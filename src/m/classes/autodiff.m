%AUTODIFF class definition
%
%   Usage:
%      autodiff=autodiff();

classdef autodiff
	properties (SetAccess=public)  
		% {{{ 
		isautodiff            = false;
		dependents            = {};
		independents          = {};
		driver                = 'fos_forward';
		obufsize              = NaN;
		lbufsize              = NaN;
		cbufsize              = NaN;
		tbufsize              = NaN;
		gcTriggerRatio        = NaN;
		gcTriggerMaxSize      = NaN;
		tapeAlloc             = NaN;
		outputTapeMemory      = 0;
		outputTime            = 0;
		enablePreaccumulation = 0;
		end
		%}}}
	methods
		function self = autodiff(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{
		self.obufsize         = 524288;
		self.lbufsize         = 524288;
		self.cbufsize         = 524288;
		self.tbufsize         = 524288;
		self.gcTriggerRatio   = 2.0;
		self.gcTriggerMaxSize = 65536;
		self.tapeAlloc        = 15000000;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return 
			if ~self.isautodiff, return; end

			%Driver value:
			md = checkfield(md,'fieldname','autodiff.driver','values',{'fos_forward','fov_forward','fov_forward_all','fos_reverse','fov_reverse','fov_reverse_all'});
			
			%buffer values: 
			md = checkfield(md,'fieldname','autodiff.obufsize','>=',16);
			md = checkfield(md,'fieldname','autodiff.lbufsize','>=',16);
			md = checkfield(md,'fieldname','autodiff.cbufsize','>=',16);
			md = checkfield(md,'fieldname','autodiff.tbufsize','>=',16);
			md = checkfield(md,'fieldname','autodiff.gcTriggerRatio','>=',0);
			md = checkfield(md,'fieldname','autodiff.gcTriggerMaxSize','>=',65536);
			md = checkfield(md,'fieldname','autodiff.tapeAlloc','>=',0);

			% Memory and time output
			md = checkfield(md,'fieldname','autodiff.outputTapeMemory','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','autodiff.outputTime','numel',[1],'values',[0 1]);

			% Memory reduction options
			md = checkfield(md,'fieldname','autodiff.enablePreaccumulation','>=',0);

			%go through our dependents and independents and check consistency: 
			for i=1:numel(self.dependents),
				dep=self.dependents{i};
				if isempty(dep)
					md = checkmessage(md,['md.autodiff.dependents{' num2str(i) '} is empty!']);
				else
					md=checkconsistency(dep,md,solution,analyses);
				end
			end
			for i=1:numel(self.independents),
				indep=self.independents{i};
				if isempty(indep)
					md = checkmessage(md,['md.autodiff.independents{' num2str(i) '} is empty!']);
				else
					md=checkconsistency(indep,md,i,solution,analyses,self.driver);
				end
			end

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   automatic differentiation parameters:'));
			fielddisplay(self,'isautodiff','indicates if the automatic differentiation is activated');
			fielddisplay(self,'dependents','list of dependent variables');
			fielddisplay(self,'independents','list of independent variables');
			fielddisplay(self,'driver','ADOLC driver (''fos_forward'' or ''fov_forward'')');
			fielddisplay(self,'obufsize','Number of operations per buffer (==OBUFSIZE in usrparms.h)');
			fielddisplay(self,'lbufsize','Number of locations per buffer (==LBUFSIZE in usrparms.h)');
			fielddisplay(self,'cbufsize','Number of values per buffer (==CBUFSIZE in usrparms.h)');
			fielddisplay(self,'tbufsize','Number of taylors per buffer (<=TBUFSIZE in usrparms.h)');
			fielddisplay(self,'gcTriggerRatio','free location block sorting/consolidation triggered if the ratio between allocated and used locations exceeds gcTriggerRatio');
			fielddisplay(self,'gcTriggerMaxSize','free location block sorting/consolidation triggered if the allocated locations exceed gcTriggerMaxSize');
			fielddisplay(self,'tapeAlloc','Iteration count of a priori memory allocation of the AD tape');
			fielddisplay(self,'outputTapeMemory','Write AD tape memory statistics to file ad_mem.dat');
			fielddisplay(self,'outputTime','Write AD recording and evaluation times to file ad_time.dat');
			fielddisplay(self,'enablePreaccumulation','Enable CoDiPack preaccumulation in augmented places');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			WriteData(fid,prefix,'object',self,'fieldname','isautodiff','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','driver','format','String');

			%early return
			if ~self.isautodiff
				WriteData(fid,prefix,'data',false,'name','md.autodiff.mass_flux_segments_present','format','Boolean');
				WriteData(fid,prefix,'data',false,'name','md.autodiff.keep','format','Boolean');
				return;
			end

			%buffer sizes
			WriteData(fid,prefix,'object',self,'fieldname','obufsize','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','lbufsize','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','cbufsize','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','tbufsize','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','gcTriggerRatio','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','gcTriggerMaxSize','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','tapeAlloc','format','Integer');

			%output of memory and time
			WriteData(fid,prefix,'object',self,'fieldname','outputTapeMemory','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','outputTime','format','Boolean');

			%memory reduction options
			WriteData(fid,prefix,'object',self,'fieldname','enablePreaccumulation','format','Boolean');
			%process dependent variables
			num_dependent_objects=numel(self.dependents);
			WriteData(fid,prefix,'data',num_dependent_objects,'name','md.autodiff.num_dependent_objects','format','Integer');
			if(num_dependent_objects)
				names={};
				for i=1:num_dependent_objects,
					dep=self.dependents{i};
					names{i}=dep.name;
				end
				WriteData(fid,prefix,'data',names,'name','md.autodiff.dependent_object_names','format','StringArray');
			end

			%process independent variables
			num_independent_objects=numel(self.independents);
			WriteData(fid,prefix,'data',num_independent_objects,'name','md.autodiff.num_independent_objects','format','Integer');
			for i=1:num_independent_objects
				indep=self.independents{i};
				WriteData(fid,prefix,'data',indep.name,'name','md.autodiff.independent_name','format','String');
				WriteData(fid,prefix,'data',indep.min_parameters,'name','md.autodiff.independent_min_parameters','format','DoubleMat','mattype',3);
				WriteData(fid,prefix,'data',indep.max_parameters,'name','md.autodiff.independent_max_parameters','format','DoubleMat','mattype',3);
				WriteData(fid,prefix,'data',indep.control_scaling_factor,'name','md.autodiff.independent_scaling_factor','format','Double');
				WriteData(fid,prefix,'data',indep.control_size,'name','md.autodiff.independent_control_size','format','Integer');
			end

			%if driver is fos_forward, build inde
			if strcmpi(self.driver,'fos_forward')
				index=0;
				for i=1:num_independent_objects,
					indep=self.independents{i};
					if ~isnan(indep.fos_forward_index),
						index=index+indep.fos_forward_index;
						break;
					else
						if strcmpi(indep.type,'scalar'),
							index=index+1;
						else
							index=index+indep.nods;
						end
					end
				end
				index=index-1; %get c-index numbering going
				WriteData(fid,prefix,'data',index,'name','md.autodiff.fos_forward_index','format','Integer');
			end

			%if driver is fos_reverse, build index:
			if strcmpi(self.driver,'fos_reverse'),
				index=0;

				for i=1:num_dependent_objects,
					dep=self.dependents{i};
					if ~isnan(dep.fos_reverse_index),
						index=index+dep.fos_reverse_index;
						break;
					else
						index=index+1;
					end
				end
				index=index-1; %get c-index numbering going
				WriteData(fid,prefix,'data',index,'name','md.autodiff.fos_reverse_index','format','Integer');
			end

			%if driver is fov_forward, build indices
			if strcmpi(self.driver,'fov_forward'),
				indices=0;

				for i=1:num_independent_objects,
					indep=self.independents{i};
					if ~isempty(indep.fos_forward_index),
						indices=indices+indep.fov_forward_indices;
						break;
					else
						if strcmpi(indep.type,'scalar'),
							indices=indices+1;
						else
							indices=indices+indep.nods;
						end
					end
				end
				indices=indices-1; %get c-indices numbering going
				WriteData(fid,prefix,'data',indices,'name','md.autodiff.fov_forward_indices','format','IntMat','mattype',3);
			end

			%deal with mass fluxes
			mass_flux_segments=cell(0,1);
			for i=1:num_dependent_objects,
				dep=self.dependents{i};
				if strcmpi(dep.name,'MassFlux'),
					mass_flux_segments{end+1,1}=dep.segments;
				end
			end
			if ~isempty(mass_flux_segments), 
				WriteData(fid,prefix,'data',mass_flux_segments,'name','md.autodiff.mass_flux_segments','format','MatArray');
				flag=true;
			else
				flag=false;
			end
			WriteData(fid,prefix,'data',flag,'name','md.autodiff.mass_flux_segments_present','format','Boolean');

			%deal with trace keep on
			keep=false;

			%From ADOLC userdoc: 
			% The optional integer argument keep of trace on determines whether the numerical values of all active variables are 
			% recorded in a buffered temporary array or file called the taylor stack. This option takes effect if keep = 1 and 
			% prepares the scene for an immediately following gradient evaluation by a call to a routine implementing the reverse 
			% mode as described in the Section 4 and Section 5. 
			%

			if length(self.driver)<=3,
				keep=false; %there is no "_reverse" string within the driver string: 
			else
				if strncmpi(self.driver(4:end),'_reverse',8),
					keep=true;
				else
					keep=false;
				end
			end
			WriteData(fid,prefix,'data',keep,'name','md.autodiff.keep','format','Boolean');

		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
			%do nothing for now
			if self.isautodiff,
				error('autodiff savemodeljs error message: not implemented yet!');
			end
		end % }}}
	end
end
