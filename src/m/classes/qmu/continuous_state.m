%
%  definition for the continuous_state class.
%
%  [csv]=continuous_state(varargin)
%
%  where the required varargin are:
%    descriptor    (char, description, '')
%    initst        (double, initial state, 0.)
%  and the optional varargin and defaults are:
%    lower         (double, lower bound, -Inf)
%    upper         (double, upper bound,  Inf)
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
classdef continuous_state
    properties
        descriptor='';
        initst    = 0.;
        lower     =-Inf;
        upper     = Inf;
    end

    methods
        function [csv]=continuous_state(varargin)

            switch nargin

%  create a default object

                case 0

%  copy the object

                case 1
                    if isa(varargin{1},'continuous_state')
                        csv=varargin{1};
                    else
                        error('Object ''%s'' is a ''%s'' class object, not ''%s''.',...
                            inputname(1),class(varargin{1}),'continuous_state');
                    end

%  create the object from the input

                otherwise
                    asizec=num2cell(array_size(varargin{1:min(nargin,4)}));
                    csv(asizec{:})=continuous_state;
                    clear asizec

                    if ischar(varargin{1})
                        varargin{1}=cellstr(varargin{1});
                    end
                    for i=1:numel(csv)
                        if (numel(varargin{1}) > 1)
                            csv(i).descriptor=varargin{1}{i};
                        else
                            csv(i).descriptor=[char(varargin{1}) string_dim(csv,i,'vector')];
                        end
                    end

                    if (nargin >= 2)
                        for i=1:numel(csv)
                            if (numel(varargin{2}) > 1)
                                csv(i).initst    =varargin{2}(i);
                            else
                                csv(i).initst    =varargin{2};
                            end
                        end
                        if (nargin >= 3)
                            for i=1:numel(csv)
                                if (numel(varargin{3}) > 1)
                                    csv(i).lower     =varargin{3}(i);
                                else
                                    csv(i).lower     =varargin{3};
                                end
                            end
                            if (nargin >= 4)
                                for i=1:numel(csv)
                                    if (numel(varargin{4}) > 1)
                                        csv(i).upper     =varargin{4}(i);
                                    else
                                        csv(i).upper     =varargin{4};
                                    end
                                end
                                if (nargin > 4)
                                    warning('continuous_state:extra_arg',...
                                        'Extra arguments for object of class ''%s''.',...
                                        class(csv));
                                end
                            end
                        end
                    end
            end

        end

        function []=disp(csv)

%  display the object

            disp(sprintf('\n'));
            for i=1:numel(csv)
                disp(sprintf('class ''%s'' object ''%s%s'' = \n',...
                    class(csv),inputname(1),string_dim(csv,i)));
                disp(sprintf('    descriptor: ''%s'''  ,csv(i).descriptor));
                disp(sprintf('        initst: %g'      ,csv(i).initst));
                disp(sprintf('         lower: %g'      ,csv(i).lower));
                disp(sprintf('         upper: %g\n'    ,csv(i).upper));
            end

        end

        function [desc]  =prop_desc(csv,dstr)
            desc=cell(1,numel(csv));
            for i=1:numel(csv)
                if ~isempty(csv(i).descriptor)
                    desc(i)=cellstr(csv(i).descriptor);
                elseif ~isempty(inputname(1))
                    desc(i)=cellstr([inputname(1) string_dim(csv,i,'vector')]);
                elseif exist('dstr','var')
                    desc(i)=cellstr([dstr         string_dim(csv,i,'vector')]);
                else
                    desc(i)=cellstr(['csv'        string_dim(csv,i,'vector')]);
                end
            end
            desc=allempty(desc);
        end
        function [initpt]=prop_initpt(csv)
            initpt=[];
        end
        function [lower] =prop_lower(csv)
            lower=zeros(1,numel(csv));
            for i=1:numel(csv)
                lower(i)=csv(i).lower;
            end
            lower=allequal(lower,-Inf);
        end
        function [upper] =prop_upper(csv)
            upper=zeros(1,numel(csv));
            for i=1:numel(csv)
                upper(i)=csv(i).upper;
            end
            upper=allequal(upper, Inf);
        end
        function [mean]  =prop_mean(csv)
            mean=[];
        end
        function [stddev]=prop_stddev(csv)
            stddev=[];
        end
        function [initst]=prop_initst(csv)
            initst=zeros(1,numel(csv));
            for i=1:numel(csv)
                initst(i)=csv(i).initst;
            end
            initst=allequal(initst,0.);
        end
        function [stype] =prop_stype(csv)
            stype={};
        end
        function [scale] =prop_scale(csv)
            scale=[];
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


    end

    methods (Static)
        function []=dakota_write(fidi,dvar)

%  collect only the variables of the appropriate class

            csv=struc_class(dvar,'continuous_state');

%  write variables

            vlist_write(fidi,'continuous_state','csv',csv);
        end
    end
end
