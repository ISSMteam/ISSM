%NORMAL_UNCERTAIN class definition
%
%   Usage:
%      nuv=normal_uncertain('descriptor',descriptor,'mean',mean,'stddev',stddev,'partition',partition);
%      where nuv is the normal_uncertain object returned by the constructor, mean and stddev are self
%      explanatory.  partition is the partition vector for distributed variables. Can be a partition
%      vector over elements or vertices.
%
%   Example:
%      md.qmu.variables.rheology=normal_uncertain('descriptor','RheologyBBar','mean',1,'stddev',.05);
%      md.qmu.variables.rheology=normal_uncertain('descriptor','scaled_RheologyBBar','mean',1,'stddev',.05,'partition',vpartition);
%

classdef normal_uncertain
	properties
		descriptor	= '';
		mean		= NaN;
		stddev		= NaN;
		partition	= [];
		nsteps		= 0;
	end
	methods
		function self=normal_uncertain(varargin) %constructor {{{

			%recover options:
			options=pairoptions(varargin{:});

			%initialize fields:
			self.descriptor=getfieldvalue(options,'descriptor');
			self.mean=getfieldvalue(options,'mean');
			self.stddev=getfieldvalue(options,'stddev');

			%if the variable is scaled,  a partition vector should have been 
			%supplied, and that partition vector should have as many partitions 
			%as the mean and stddev vectors:
			if self.isscaled(),
				self.partition=getfieldvalue(options,'partition');
				self.nsteps=getfieldvalue(options,'nsteps',1);
				npart=qmupart2npart(self.partition);
				if npart~=size(self.mean,1),
					error(['normal_uncertain constructor: for the scaled variable ' self.descriptor ' the row size of the mean field should be identical to the number of partitions']);
				end
				if npart~=size(self.stddev,1),
					error(['normal_uncertain constructor: for the scaled variable ' self.descriptor ' the row size of the stddev field should be identical to the number of partitions']);
				end
				if self.nsteps~=size(self.mean,2),
					error(['normal_uncertain constructor: for the scaled variable ' self.descriptor ' the col size of the mean field should be identical to the number of time steps']);
				end
				if self.nsteps~=size(self.stddev,2),
					error(['normal_uncertain constructor: for the scaled variable ' self.descriptor ' the col size of the stddev field should be identical to the number of time steps']);
				end

			end

		end %}}}
		function disp(self) % {{{
			disp(sprintf('   normal uncertain variable: '));
			fielddisplay(self,'descriptor','name tag');
			fielddisplay(self,'mean','pdf mean');
			fielddisplay(self,'stddev','pdf standard deviation');
			if ~isempty(self.partition),
				fielddisplay(self,'partition','partition vector defining where sampling will occur');
			end
			fielddisplay(self,'nsteps','number of time steps');
		end
		%}}}
		function md=checkconsistency(self,md,solution,analyses) % {{{

			md = checkfield(md,'field',self.mean,'fieldname','normal_uncertain.mean','NaN',1,'Inf',1,'>=',0);
			md = checkfield(md,'field',self.stddev,'fieldname','normal_uncertain.stddev','NaN',1,'Inf',1,'>=',0);
			if self.isscaled(),
				if isempty(self.partition),
					error('normal_uncertain is a scaled variable, but it''s missing a partition vector');
				end
				%better have a partition vector that has as many partitions as stddev's size:
				if size(self.stddev,1)~=partition_npart(self.partition),
					error('normal_uncertain error message: row size of stddev and partition size should be identical');
				end
				if size(self.mean,1)~=partition_npart(self.partition),
					error('normal_uncertain error message: row size of mean and partition size should be identical');
				end
				%we need as steps in stddev and mean as there are time steps: 
				if size(self.stddev,2)~=self.nsteps,
					error('normal_uncertain error message: col size of stddev and number of time steps should be identical');
				end
				if size(self.mean,2)~=self.nsteps,
					error('normal_uncertain error message: col size of mean and number of time steps should be identical');
				end

				md = checkfield(md,'field',self.partition,'fieldname','normal_uncertain.partition','NaN',1,'Inf',1,'>=',-1,'numel',[md.mesh.numberofvertices,md.mesh.numberofelements]);
				if size(self.partition,2)>1,
					error('normal_uncertain error message: partition should be a column vector');
				end
				partcheck=unique(self.partition);
				partmin=min(partcheck);
				partmax=max(partcheck);
				if partmax<-1,
					error('normal_uncertain error message: partition vector''s min value should be -1 (for no partition), or start at 0');
				end
				nmax=max(md.mesh.numberofelements,md.mesh.numberofvertices);
				if partmax>nmax,
					error('normal_uncertain error message: partition vector''s values cannot go over the number of vertices or elements');
				end
			end
		end % }}}
		%virtual functions needed by qmu processing algorithms:
		%implemented:
		function [desc]=prop_desc(nuv,dstr) % {{{
			desc=cell(1,numel(nuv));
			for i=1:numel(nuv)
				if ~isempty(nuv(i).descriptor)
					desc(i)=cellstr(nuv(i).descriptor);
				elseif ~isempty(inputname(1))
					desc(i)=cellstr([inputname(1) string_dim(nuv,i,'vector')]);
				elseif exist('dstr','var')
					desc(i)=cellstr([dstr         string_dim(nuv,i,'vector')]);
				else
					desc(i)=cellstr(['nuv'        string_dim(nuv,i,'vector')]);
				end
			end
			desc=allempty(desc);
		end %}}}
		function [mean]=prop_mean(nuv) % {{{
			mean=zeros(1,numel(nuv));
			for i=1:numel(nuv)
				mean(i)=nuv(i).mean;
			end
		end % }}}
		function [stddev]=prop_stddev(nuv) % {{{
			stddev=zeros(1,numel(nuv));
			for i=1:numel(nuv)
				stddev(i)=nuv(i).stddev;
			end
		end % }}}
		function [lower]=prop_lower(nuv) % {{{
			lower=[];
		end % }}}
		function [upper]=prop_upper(nuv) % {{{
			upper=[];
		end % }}}
		function [abscissas]=prop_abscissas(hbu) % {{{
			abscissas=[];
		end % }}}
		function [counts]=prop_counts(hbu) % {{{
			counts=[];
		end % }}}
		function [pairs_per_variable]=prop_pairs_per_variable(hbu) % {{{
			pairs_per_variable=[];
		end % }}}
		function [initpt]=prop_initpt(nuv) % {{{
			initpt=[];
		end % }}}
		function [initst]=prop_initst(nuv) % {{{
			initst=[];
		end % }}}
		function [stype]=prop_stype(nuv) % {{{
			stype={};
		end % }}}
		function [scale]=prop_scale(nuv) % {{{
			scale=[];
		end % }}}
		%new methods:
		function distributed=isdistributed(self) % {{{
			if strncmp(self.descriptor,'distributed_',12),
				distributed=1;
			else
				distributed=0;
			end
		end % }}}
		function scaled=isscaled(self) % {{{
			if strncmp(self.descriptor,'scaled_',7),
				scaled=1;
			else
				scaled=0;
			end
		end % }}}
	end
	methods (Static)
		function []=dakota_write(fidi,dvar) % {{{
			%  collect only the variables of the appropriate class
			nuv=struc_class(dvar,'normal_uncertain');
			%  write variables
			vlist_write(fidi,'normal_uncertain','nuv',nuv);
		end % }}}
	end
end
