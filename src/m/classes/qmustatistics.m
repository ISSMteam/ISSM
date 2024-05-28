%QMUSTATISTICS class definition
%
%   Usage:
%      stats=qmustatistics();

classdef qmustatistics
	properties (SetAccess=public)
	
		%generic information: 
		nfiles_per_directory = 5;  %number of files per output directory.
		ndirectories         = 50; %number of output directories; should be < numcpus
	
		method = struct('name',{'None'},'fields',{[]},'steps',{[]},'nbins',{NaN},'indices',{[]});

		%name: name of method, one of 'None','Histogram','SampleSeries' or 'MeanVariance'
		%fields: fields for the  statistics being requested, ex: 'Sealevel','BslcIce','BslcHydro'
		%steps: time steps at which each field statistic is computed, ex: [1:2,5,20] or 1:100
		%nbins: number of bins for 'Histogram' statistics
		%indices: vertex indices at which to retrieve samples

	end
	methods
		function self = qmustatistics(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					inputstruct=varargin{1};
					list1 = properties('qmustatistics');
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

			self.method.name='None';
			self.nfiles_per_directory = 5;  %number of files per output directory.
			self.ndirectories         = 50; %number of output directories.  %ndirectories should be < numcpus

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return: 
			if md.qmu.isdakota==0; return; end;
			if strcmpi(self.method(1).name,'None'),return;end;

			%Checks:
			md = checkfield(md,'fieldname','qmu.statistics.nfiles_per_directory','>=',1);
			if self.ndirectories>md.cluster.np, 
				error('qmustatistics consistency check: number of cluster CPUs should be > number of output directories');
			end
			if (self.ndirectories*self.nfiles_per_directory)~=md.qmu.method.params.samples,
				error('qmustatistics consistency check: number of directories x number of files per directory should be == to number of samples requested!');
			end
			for i=1:length(self.method),
				m=self.method(i);
				if strcmpi(m.name,'Histogram'),
					md = checkfield(md,'fieldname',sprintf('qmu.statistics.method(%i).nbins',i),'>=',1,'<=',md.qmu.method.params.samples);
				end
				for f=1:length(m.fields),
					if ~ischar(m.fields{f}),
						error(sprintf('qmustatistics consistency check error: qmu.statistics.method(%i).fields{%i} is not a string!',i,f));
					end
				end
				for s=1:length(m.steps),
					if m.steps(s)<=0,
						error(sprintf('qmustatistics consistency check error: qmu.statistics.method(%i).steps(%i) should be > 0!',i,s));
					end
					if m.steps(s)> md.mesh.numberofvertices
						error(sprintf('qmustatistics consistency check error: qmu.statistics.method(%i).steps(%i) should be < md.mesh.numberofvertices!',i,s));
					end
				end

			end

		end % }}}
		function list=defaultoutputs(self,md) % {{{
		end % }}}
		function disp(self) % {{{

			disp(sprintf('qmustatistics: post-Dakota run processing of QMU statistics:'));
				
			if strcmpi(self.method(1).name,'None'),  return; end;

			%generic information:
			fielddisplay(self,'nfiles_per_directory','number of files per output directory');
			fielddisplay(self,'ndirectories','number of output directories; should be < numcpus');

			for i=1:length(self.method),
				disp(['   method #' num2str(i)]);
				disp(self.method(i));
			end

		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			if strcmpi(self.method(1).name,'None'), 
				WriteData(fid,prefix,'name','md.qmu.statistics','data',0,'format','Boolean');
				statistics=0; 
				return;
			else 
				WriteData(fid,prefix,'name','md.qmu.statistics','data',1,'format','Boolean');
				statistics=1;
			end

			if statistics,
				WriteData(fid,prefix,'name','md.qmu.statistics.nfiles_per_directory','data',self.nfiles_per_directory,'format','Integer');
				WriteData(fid,prefix,'name','md.qmu.statistics.ndirectories','data',self.ndirectories,'format','Integer');
				WriteData(fid,prefix,'name','md.qmu.statistics.numstatistics','data',length(self.method),'format','Integer');
				for i=1:length(self.method),
					m=self.method(i); 
					WriteData(fid,prefix,'name',sprintf('md.qmu.statistics.method(%i).name',i),'data',m.name,'format','String');
					WriteData(fid,prefix,'data',m.fields,'name',sprintf('md.qmu.statistics.method(%i).fields',i),'format','StringArray');
					WriteData(fid,prefix,'data',m.steps,'name',sprintf('md.qmu.statistics.method(%i).steps',i),'format','IntMat','mattype',3);

					if strcmpi(m.name,'Histogram'),
						WriteData(fid,prefix,'name',sprintf('md.qmu.statistics.method(%i).nbins',i),'data',m.nbins,'format','Integer');
					elseif strcmpi(m.name,'MeanVariance'),
						%do nothing
					elseif strcmpi(m.name,'SampleSeries'),
						WriteData(fid,prefix,'data',m.indices,'name',sprintf('md.qmu.statistics.method(%i).indices',i),'format','IntMat','mattype',3);
					else 
						error(sprintf('qmustatistics marshall error message: unknown type ''%s'' for qmu.statistics.method(%i)',m.name,i));
					end

				end
			end

		end % }}}
		function self = extrude(self,md) % {{{

		end % }}}
		function savemodeljs(self,fid,modelname) % {{{

			writejs1Darray(fid,[modelname '.stressbalance.spcvx'],self.spcvx);
			writejs1Darray(fid,[modelname '.stressbalance.spcvy'],self.spcvy);
			writejs1Darray(fid,[modelname '.stressbalance.spcvz'],self.spcvz);
			writejsdouble(fid,[modelname '.stressbalance.restol'],self.restol);
			writejsdouble(fid,[modelname '.stressbalance.reltol'],self.reltol);
			writejsdouble(fid,[modelname '.stressbalance.abstol'],self.abstol);
			writejsdouble(fid,[modelname '.stressbalance.isnewton'],self.isnewton);
			writejsdouble(fid,[modelname '.stressbalance.FSreconditioning'],self.FSreconditioning);
			writejsdouble(fid,[modelname '.stressbalance.maxiter'],self.maxiter);
			writejsdouble(fid,[modelname '.stressbalance.shelf_dampening'],self.shelf_dampening);
			writejs1Darray(fid,[modelname '.stressbalance.vertex_pairing'],self.vertex_pairing);
			writejsdouble(fid,[modelname '.stressbalance.penalty_factor'],self.penalty_factor);
			writejsdouble(fid,[modelname '.stressbalance.rift_penalty_lock'],self.rift_penalty_lock);
			writejsdouble(fid,[modelname '.stressbalance.rift_penalty_threshold'],self.rift_penalty_threshold);
			writejs2Darray(fid,[modelname '.stressbalance.referential'],self.referential);
			writejs2Darray(fid,[modelname '.stressbalance.loadingforce'],self.loadingforce);
			writejscellstring(fid,[modelname '.stressbalance.requested_outputs'],self.requested_outputs);

		end % }}}
	end
end
