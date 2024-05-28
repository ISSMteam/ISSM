%UNIFORM_UNCERTAIN Class definition
%
%%   Usage:
%      nuv=uniform_uncertain('descriptor',descriptor,'lower',lower,'upper',upper,'partition',partition);
%      where nuv is the uniform_uncertain object returned by the constructor, lower and upper are the 
%      pdf distribution bounds, and partition is the partition vector for distributed variables. 
%      Can be a partition %      vector over elements or vertices.
% 
%   Example:
%      md.qmu.variables.rheology=uniform_uncertain('descriptor','RheologyBBar','lower',1e8,'upper',1e9);
%      md.qmu.variables.rheology=uniform_uncertain('descriptor','RheologyBBar','lower',1e8,'upper',1e9,'partition',vpartition);
% 

classdef uniform_uncertain
	properties
		descriptor	= '';
		lower		= -Inf;
		upper		= Inf;
		partition	= [];
		nsteps		= 0;
	end
	methods
		function self=uniform_uncertain(varargin) %constructor {{{

			%recover options: 
			options=pairoptions(varargin{:});

			%initialize fields: 
			self.descriptor=getfieldvalue(options,'descriptor');
			self.upper=getfieldvalue(options,'upper');
			self.lower=getfieldvalue(options,'lower');

			%if the variable is scaled, a partition vector should have been 
			%supplied, and that partition vector should have as many partitions 
			%as the lower and upper vectors:
			if self.isscaled(),
				self.partition=getfieldvalue(options,'partition');
				self.nsteps=getfieldvalue(options,'nsteps',1);
				npart=qmupart2npart(self.partition);
				if npart~=size(self.upper,1),
					error(['uniform_uncertain constructor: for the scaled variable' self.descriptor ' the row size of the upper field should be identical to the number of partitions']);
				end
				if npart~=size(self.lower,1),
					error(['uniform_uncertain constructor: for the scaled variable' self.descriptor ' the row size of the lower field should be identical to the number of partitions']);
				end
				if self.nsteps~=size(self.upper,2),
					error(['uniform_uncertain constructor: for the scaled variable ' self.descriptor ' the col size of the upper field should be identical to the number of time steps']);
				end
				if self.nsteps~=size(self.lower,2),
					error(['uniform_uncertain constructor: for the scaled variable ' self.descriptor ' the col size of the lower field should be identical to the number of time steps']);
				end
			end

		end %}}}
		function disp(self) % {{{

			disp(sprintf('   uniform uncertain variable: '));
			fielddisplay(self,'descriptor','name tag');
			fielddisplay(self,'lower','pdf lower bound');
			fielddisplay(self,'upper','pdf upper bound');
			if ~isempty(self.partition),
				fielddisplay(self,'partition','partition vector defining where sampling will occur');
			end
			fielddisplay(self,'nsteps','number of time steps');
		end 
		%}}}
		function md=checkconsistency(self,md,solution,analyses) % {{{

			md = checkfield(md,'field',self.upper,'fieldname','uniform_uncertain.upper','NaN',1,'Inf',1,'>',self.lower,'numel',length(self.lower));
			md = checkfield(md,'field',self.lower,'fieldname','uniform_uncertain.lower','NaN',1,'Inf',1,'<',self.upper,'numel',length(self.upper));
			if self.isscaled(),
				if isempty(self.partition),
					error('uniform_uncertain is a scaled variable, but it''s missing a partition vector');
				end
				%better have a partition vector that has as many partitions as upper and lower's size: 
				if size(self.upper,1)~=partition_npart(self.partition),
					error('uniform_uncertain error message: row size of upper and partition size should be identical');
				end
				if size(self.lower,1)~=partition_npart(self.partition),
					error('uniform_uncertain error message: row size of lower and partition size should be identical');
				end
				%we need as steps in upper and lower as there are time steps: 
				if size(self.upper,2)~=self.nsteps,
					error('uniform_uncertain error message: col size of upper and number of time steps should be identical');
				end
				if size(self.lower,2)~=self.nsteps,
					error('uniform_uncertain error message: col size of lower and number of time steps should be identical');
				end

				md = checkfield(md,'field',self.partition,'fieldname','uniform_uncertain.partition','NaN',1,'Inf',1,'>=',-1,'numel',[md.mesh.numberofvertices,md.mesh.numberofelements]);
				if size(self.partition,2)>1,
					error('uniform_uncertain error message: partition should be a column vector');
				end
				partcheck=unique(self.partition); 
				partmin=min(partcheck); 
				partmax=max(partcheck);
				if partmax<-1,
					error('uniform_uncertain error message: partition vector''s min value should be -1 (for no partition), or start at 0');
				end
				nmax=max(md.mesh.numberofelements,md.mesh.numberofvertices);
				if partmax>nmax,
					error('uniform_uncertain error message: partition vector''s values cannot go over the number of vertices or elements');
				end
			end
		end % }}}
		%virtual functions needed by qmu processing algorithms:
		%implemented:
		function [desc]=prop_desc(uuv,dstr) % {{{
			desc=cell(1,numel(uuv));
			for i=1:numel(uuv)
				if ~isempty(uuv(i).descriptor)
					desc(i)=cellstr(uuv(i).descriptor);
				elseif ~isempty(inputname(1))
					desc(i)=cellstr([inputname(1) string_dim(uuv,i,'vector')]);
				elseif exist('dstr','var')
					desc(i)=cellstr([dstr         string_dim(uuv,i,'vector')]);
				else
					desc(i)=cellstr(['uuv'        string_dim(uuv,i,'vector')]);
				end
			end
			desc=allempty(desc);
		end %}}}
		function [lower]=prop_lower(uuv) % {{{
			lower=zeros(1,numel(uuv));
			for i=1:numel(uuv)
				lower(i)=uuv(i).lower;
			end
			lower=allequal(lower,-Inf);
		end %}}}
		function [upper]=prop_upper(uuv) % {{{
			upper=zeros(1,numel(uuv));
			for i=1:numel(uuv)
				upper(i)=uuv(i).upper;
			end
			%upper=allequal(upper, Inf);
		end % }}}
		%default
		function [stddev]=prop_stddev(uuv)  %{{{
			stddev=[];
		end % }}}
		function [mean]=prop_mean(nuv) % {{{
			mean=[];
		end % }}}
		function [initpt]=prop_initpt(uuv) %{{{
			initpt=[];
		end %}}}
		function [initst]=prop_initst(uuv) %{{{
			initst=[];
		end %}}}
		function [stype]=prop_stype(uuv) %{{{
			stype={};
		end %}}}
		function [scale]=prop_scale(uuv) %{{{
			scale=[];
		end %}}}
		function [abscissas]=prop_abscissas(hbu) % {{{
			abscissas=[]; 
		end % }}}
		function [counts]=prop_counts(hbu) % {{{
			counts=[]; 
		end % }}}
		function [pairs_per_variable]=prop_pairs_per_variable(hbu) % {{{
			pairs_per_variable=[];
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
			uuv=struc_class(dvar,'uniform_uncertain');
			%  write variables
			vlist_write(fidi,'uniform_uncertain','uuv',uuv);
		end %}}}
	end
end
