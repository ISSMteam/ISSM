%
%  definition for the continuous_design class.
%
%  [cdv]=continuous_design(varargin)
%
%  where the required varargin are:
%    descriptor    (char, description, '')
%    initpt        (double, initial point, 0.)
%  and the optional varargin and defaults are:
%    lower         (double, lower bound, -Inf)
%    upper         (double, upper bound,  Inf)
%    scale_type    (char, scaling type, 'none')
%    scale         (double, scaling factor, 1.)
%
%  note that zero arguments constructs a default instance; one
%  argument of the class copies the instance; and two or more
%  arguments constructs a new instance from the arguments.
%
%  "Copyright 2009, by the California Institute of Technology.
%  ALL RIGHTS RESERVED. United States Government Sponsorship
%  acknowledged. Any commercial use must be negotiated with
%  the Office of Technology Transfer at the California Institute
%  of Technology.  (J. Schiermeier, NTR 47078)
%
%  This software may be subject to U.S. export control laws.
%  By accepting this  software, the user agrees to comply with
%  all applicable U.S. export laws and regulations. User has the
%  responsibility to obtain export licenses, or other export
%  authority as may be required before exporting such information
%  to foreign countries or providing access to foreign persons."
%
classdef continuous_design
	properties
		descriptor='';
		initpt    = 0.;
		lower     =-Inf;
		upper     = Inf;
		scale_type='none';
		scale     = 1.;
		partition = [];
		nsteps    = 0;
	end

	methods
		function self=continuous_design(varargin) % {{{

			%recover options:
			options=pairoptions(varargin{:});

			%initialize fields:
			self.descriptor=getfieldvalue(options,'descriptor');
			self.initpt=getfieldvalue(options,'initpt');
			self.lower=getfieldvalue(options,'lower');
			self.upper=getfieldvalue(options,'upper');
			self.scale_type=getfieldvalue(options,'scale_type','none');
			self.scale=getfieldvalue(options,'scale',1.);

			%if the variable is scaled,  a partition vector should have been 
			%supplied, and that partition vector should have as many partitions 
			%as the upper and lower vectors: 
			if self.isscaled() | self.isdistributed(),
				self.partition=getfieldvalue(options,'partition');
				self.nsteps=getfieldvalue(options,'nsteps',1);
				npart=qmupart2npart(self.partition);
				if npart~=size(self.upper,1),
					error(['continuous_design constructor: for the scaled variable ' self.descriptor ' the row size of the upper field should be identical to the number of partitions']);
				end
				if npart~=size(self.lower,1),
					error(['continuous_design constructor: for the scaled variable ' self.descriptor ' the row size of the lower field should be identical to the number of partitions']);
				end
				if self.nsteps~=size(self.upper,2),
					error(['continuous_design constructor: for the scaled variable ' self.descriptor ' the col size of the upper field should be identical to the number of time steps']);
				end
				if self.nsteps~=size(self.lower,2),
					error(['continuous_design constructor: for the scaled variable ' self.descriptor ' the col size of the lower field should be identical to the number of time steps']);
				end

			end



		end % }}}
		function disp(self) % {{{

			disp(sprintf('   continuous design variable: '));
			fielddisplay(self,'descriptor','name tag');
			fielddisplay(self,'initpt','initial point');
			fielddisplay(self,'lower','lower values');
			fielddisplay(self,'upper','upper values');
			fielddisplay(self,'scale_type','scale type');
			fielddisplay(self,'scale','scale');

		end % }}}
			function md=checkconsistency(self,md,solution,analyses) % {{{

			md = checkfield(md,'field',self.upper,'fieldname','continuous_design.upper','NaN',1,'Inf',1,'>=',0);
			md = checkfield(md,'field',self.lower,'fieldname','continuous_design.lower','NaN',1,'Inf',1,'>=',0);
			if self.isscaled(),
				if isempty(self.partition),
					error('continuous_design is a scaled variable, but it''s missing a partition vector');
				end
				%better have a partition vector that has as many partitions as loer's size:
				if size(self.lower,1)~=partition_npart(self.partition),
					error('continuous_design error message: row size of lower and partition size should be identical');
				end
				if size(self.upper,1)~=partition_npart(self.partition),
					error('continuous_design error message: row size of upper and partition size should be identical');
				end
				%we need as steps in lower and upper as there are time steps: 
				if size(self.lower,2)~=self.nsteps,
					error('continuous_design error message: col size of lower and number of time steps should be identical');
				end
				if size(self.upper,2)~=self.nsteps,
					error('continuous_design error message: col size of upper and number of time steps should be identical');
				end

				md = checkfield(md,'field',self.partition,'fieldname','continuous_design.partition','NaN',1,'Inf',1,'>=',-1,'numel',[md.mesh.numberofvertices,md.mesh.numberofelements]);
				if size(self.partition,2)>1,
					error('continuous_design error message: partition should be a column vector');
				end
				partcheck=unique(self.partition);
				partmin=min(partcheck);
				partmax=max(partcheck);
				if partmax<-1,
					error('continuous_design error message: partition vector''s min value should be -1 (for no partition), or start at 0');
				end
				nmax=max(md.mesh.numberofelements,md.mesh.numberofvertices);
				if partmax>nmax,
					error('continuous_design error message: partition vector''s values cannot go over the number of vertices or elements');
				end
			end
		end % }}}
		function [desc]  =prop_desc(cdv,dstr) % {{{
			desc=cell(1,numel(cdv));
			for i=1:numel(cdv)
				if ~isempty(cdv(i).descriptor)
					desc(i)=cellstr(cdv(i).descriptor);
				elseif ~isempty(inputname(1))
					desc(i)=cellstr([inputname(1) string_dim(cdv,i,'vector')]);
				elseif exist('dstr','var')
					desc(i)=cellstr([dstr         string_dim(cdv,i,'vector')]);
				else
					desc(i)=cellstr(['cdv'        string_dim(cdv,i,'vector')]);
				end
			end
			desc=allempty(desc);
		end % }}}
		function [initpt]=prop_initpt(cdv) % {{{
			initpt=zeros(1,numel(cdv));
			for i=1:numel(cdv)
				initpt(i)=cdv(i).initpt;
			end
			initpt=allequal(initpt,0.);
		end % }}}
		function [lower] =prop_lower(cdv) % {{{
			lower=zeros(1,numel(cdv));
			for i=1:numel(cdv)
				lower(i)=cdv(i).lower;
			end
			lower=allequal(lower,-Inf);
		end % }}}
		function [upper] =prop_upper(cdv) % {{{
			upper=zeros(1,numel(cdv));
			for i=1:numel(cdv)
				upper(i)=cdv(i).upper;
			end
			upper=allequal(upper, Inf);
		end % }}}
		function [mean]  =prop_mean(cdv) % {{{
			mean=[];
		end % }}}
		function [stddev]=prop_stddev(cdv)  % {{{
			stddev=[];
		end % }}}
		function [initst]=prop_initst(cdv) % {{{
			initst=[];
		end % }}}
		function [stype] =prop_stype(cdv) % {{{
			stype=cell(1,numel(cdv));
			for i=1:numel(cdv)
				stype(i)=cellstr(cdv(i).scale_type);
			end
			stype=allequal(stype,'none');
		end % }}}
		function [scale] =prop_scale(cdv) % {{{
			scale=zeros(1,numel(cdv));
			for i=1:numel(cdv)
				scale(i)=cdv(i).scale;
			end
			scale=allequal(scale,1.);
		end % }}}
		function [abscissas] =prop_abscissas(hbu) % {{{
            abscissas=[]; 
        end % }}}
        function [counts] =prop_counts(hbu) % {{{
            counts=[]; 
        end % }}}
        function [pairs_per_variable] =prop_pairs_per_variable(hbu) % {{{
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
			cdv=struc_class(dvar,'continuous_design');
			%  write variables
			vlist_write(fidi,'continuous_design','cdv',cdv);
		end % }}}
	end
end
