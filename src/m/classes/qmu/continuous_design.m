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
	end

	methods
		function [cdv]=continuous_design(varargin)

			switch nargin

				%  create a default object

				case 0

					%  copy the object

				case 1
					if isa(varargin{1},'continuous_design')
						cdv=varargin{1};
					else
						error('Object ''%s'' is a ''%s'' class object, not ''%s''.',...
							inputname(1),class(varargin{1}),'continuous_design');
					end

					%  create the object from the input

				otherwise
					asizec=num2cell(array_size(varargin{1:min(nargin,6)}));
					cdv(asizec{:})=continuous_design;
					clear asizec

					if ischar(varargin{1})
						varargin{1}=cellstr(varargin{1});
					end
					for i=1:numel(cdv)
						if (numel(varargin{1}) > 1)
							cdv(i).descriptor=varargin{1}{i};
						else
							cdv(i).descriptor=[char(varargin{1}) string_dim(cdv,i,'vector')];
						end
					end

					if (nargin >= 2)
						for i=1:numel(cdv)
							if (numel(varargin{2}) > 1)
								cdv(i).initpt    =varargin{2}(i);
							else
								cdv(i).initpt    =varargin{2};
							end
						end
						if (nargin >= 3)
							for i=1:numel(cdv)
								if (numel(varargin{3}) > 1)
									cdv(i).lower     =varargin{3}(i);
								else
									cdv(i).lower     =varargin{3};
								end
							end
							if (nargin >= 4)
								for i=1:numel(cdv)
									if (numel(varargin{4}) > 1)
										cdv(i).upper     =varargin{4}(i);
									else
										cdv(i).upper     =varargin{4};
									end
								end
								if (nargin >= 5)
									if ischar(varargin{5})
										varargin{5}=cellstr(varargin{5});
									end
									for i=1:numel(cdv)
										if (numel(varargin{5}) > 1)
											cdv(i).scale_type=varargin{5}{i};
										else
											cdv(i).scale_type=char(varargin{5});
										end
									end
									if (nargin >= 6)
										for i=1:numel(cdv)
											if (numel(varargin{6}) > 1)
												cdv(i).scale     =varargin{6}(i);
											else
												cdv(i).scale     =varargin{6};
											end
										end
										if (nargin > 6)
											warning('continuous_design:extra_arg',...
												'Extra arguments for object of class ''%s''.',...
												class(cdv));
										end
									end
								end
							end
						end
					end
			end

		end

		function []=disp(cdv)

			%  display the object

			disp(sprintf('\n'));
			for i=1:numel(cdv)
				disp(sprintf('class ''%s'' object ''%s%s'' = \n',...
					class(cdv),inputname(1),string_dim(cdv,i)));
				disp(sprintf('    descriptor: ''%s'''  ,cdv(i).descriptor));
				disp(sprintf('        initpt: %g'      ,cdv(i).initpt));
				disp(sprintf('         lower: %g'      ,cdv(i).lower));
				disp(sprintf('         upper: %g'      ,cdv(i).upper));
				disp(sprintf('    scale_type: ''%s'''  ,cdv(i).scale_type));
				disp(sprintf('         scale: %g\n'    ,cdv(i).scale));
			end

		end

		function [desc]  =prop_desc(cdv,dstr)
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
		end
		function [initpt]=prop_initpt(cdv)
			initpt=zeros(1,numel(cdv));
			for i=1:numel(cdv)
				initpt(i)=cdv(i).initpt;
			end
			initpt=allequal(initpt,0.);
		end
		function [lower] =prop_lower(cdv)
			lower=zeros(1,numel(cdv));
			for i=1:numel(cdv)
				lower(i)=cdv(i).lower;
			end
			lower=allequal(lower,-Inf);
		end
		function [upper] =prop_upper(cdv)
			upper=zeros(1,numel(cdv));
			for i=1:numel(cdv)
				upper(i)=cdv(i).upper;
			end
			upper=allequal(upper, Inf);
		end
		function [mean]  =prop_mean(cdv)
			mean=[];
		end
		function [stddev]=prop_stddev(cdv)
			stddev=[];
		end
		function [initst]=prop_initst(cdv)
			initst=[];
		end
		function [stype] =prop_stype(cdv)
			stype=cell(1,numel(cdv));
			for i=1:numel(cdv)
				stype(i)=cellstr(cdv(i).scale_type);
			end
			stype=allequal(stype,'none');
		end
		function [scale] =prop_scale(cdv)
			scale=zeros(1,numel(cdv));
			for i=1:numel(cdv)
				scale(i)=cdv(i).scale;
			end
			scale=allequal(scale,1.);
		end
		function [abscissas] =prop_abscissas(hbu) % {{{
            abscissas=[]; 
        end % }}}
        function [counts] =prop_counts(hbu) % {{{
            counts=[]; 
        end % }}}
        function [pairs_per_variable] =prop_pairs_per_variable(hbu) % {{{
			pairs_per_variable=[];
        end % }}}
		function checkconsistency(self,md,solution,analyses) % {{{
		end % }}}

	end

	methods (Static)
		function []=dakota_write(fidi,dvar)

			%  collect only the variables of the appropriate class

			cdv=struc_class(dvar,'continuous_design');

			%  write variables

			vlist_write(fidi,'continuous_design','cdv',cdv);
		end
	end
end
