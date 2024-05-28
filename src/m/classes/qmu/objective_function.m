%
%  definition for the objective_function class.
%
%  [of]=objective_function(varargin)
%
%  where the required varargin are:
%    descriptor    (char, description, '')
%  and the optional varargin and defaults are:
%    scale_type    (char, scaling type, 'none')
%    scale         (double, scaling factor, 1.)
%    weight        (double, weighting factor, 1.)
%
%  note that zero arguments constructs a default instance; one
%  argument of the class copies the instance; and one or more
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
classdef objective_function
    properties
        descriptor='';
        scale_type='none';
        scale     = 1.;
        weight    = 1.;
    end

    methods
        function [of]=objective_function(varargin)

            switch nargin

%  create a default object

                case 0

%  copy the object or create the object from the input

                otherwise
                    if  (nargin == 1) && isa(varargin{1},'objective_function')
                        of=varargin{1};
                    else
                        asizec=num2cell(array_size(varargin{1:min(nargin,4)}));
                        of(asizec{:})=objective_function;
                        clear asizec

                        if ischar(varargin{1})
                            varargin{1}=cellstr(varargin{1});
                        end
                        for i=1:numel(of)
                            if (numel(varargin{1}) > 1)
                                of(i).descriptor=varargin{1}{i};
                            else
                                of(i).descriptor=[char(varargin{1}) string_dim(of,i,'vector')];
                            end
                        end

                        if (nargin >= 2)
                            if ischar(varargin{2})
                                varargin{2}=cellstr(varargin{2});
                            end
                            for i=1:numel(of)
                                if (numel(varargin{2}) > 1)
                                    of(i).scale_type=varargin{2}{i};
                                else
                                    of(i).scale_type=char(varargin{2});
                                end
                            end
                            if (nargin >= 3)
                                for i=1:numel(of)
                                    if (numel(varargin{3}) > 1)
                                        of(i).scale     =varargin{3}(i);
                                    else
                                        of(i).scale     =varargin{3};
                                    end
                                end
                                if (nargin >= 4)
                                    for i=1:numel(of)
                                        if (numel(varargin{4}) > 1)
                                            of(i).weight    =varargin{4}(i);
                                        else
                                            of(i).weight    =varargin{4};
                                        end
                                    end

                                    if (nargin > 4)
                                        warning('objective_function:extra_arg',...
                                            'Extra arguments for object of class ''%s''.',...
                                            class(of));
                                    end
                                end
                            end
                        end
                    end
            end

        end

        function []=disp(of)

%  display the object

            disp(sprintf('\n'));
            for i=1:numel(of)
                disp(sprintf('class ''%s'' object ''%s%s'' = \n',...
                    class(of),inputname(1),string_dim(of,i)));
                disp(sprintf('    descriptor: ''%s'''  ,of(i).descriptor));
                disp(sprintf('    scale_type: ''%s'''  ,of(i).scale_type));
                disp(sprintf('         scale: %g'      ,of(i).scale));
                disp(sprintf('        weight: %g\n'    ,of(i).weight));
            end

        end

        function [desc]  =prop_desc(of,dstr)
            desc=cell(1,numel(of));
            for i=1:numel(of)
                if ~isempty(of(i).descriptor)
                    desc(i)=cellstr(of(i).descriptor);
                elseif ~isempty(inputname(1))
                    desc(i)=cellstr([inputname(1) string_dim(of,i,'vector')]);
                elseif exist('dstr','var')
                    desc(i)=cellstr([dstr         string_dim(of,i,'vector')]);
                else
                    desc(i)=cellstr(['of'         string_dim(of,i,'vector')]);
                end
            end
            desc=allempty(desc);
        end
        function [stype] =prop_stype(of)
            stype=cell(1,numel(of));
            for i=1:numel(of)
                stype(i)=cellstr(of(i).scale_type);
            end
            stype=allequal(stype,'none');
        end
        function [scale] =prop_scale(of)
            scale=zeros(1,numel(of));
            for i=1:numel(of)
                scale(i)=of(i).scale;
            end
            scale=allequal(scale,1.);
        end
        function [weight]=prop_weight(of)
            weight=zeros(1,numel(of));
            for i=1:numel(of)
                weight(i)=of(i).weight;
            end
            weight=allequal(weight,1.);
        end
        function [lower] =prop_lower(of)
            lower=[];
        end
        function [upper] =prop_upper(of)
            upper=[];
        end
        function [target]=prop_target(of)
            target=[];
        end
    end

    methods (Static)
        function [rdesc]=dakota_write(fidi,dresp,rdesc)

%  collect only the responses of the appropriate class

            of=struc_class(dresp,'objective_function');

%  write responses

            [rdesc]=rlist_write(fidi,'objective_functions','objective_function',of,rdesc);
        end

        function []=dakota_rlev_write(fidi,dresp,params)
        end
    end
end
