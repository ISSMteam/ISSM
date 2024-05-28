%
%  definition for the least_squares_term class.
%
%  [lst]=least_squares_term(varargin)
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
classdef least_squares_term
    properties
        descriptor='';
        scale_type='none';
        scale     = 1.;
        weight    = 1.;
    end

    methods
        function [lst]=least_squares_term(varargin)

            switch nargin

%  create a default object

                case 0

%  copy the object or create the object from the input

                otherwise
                    if  (nargin == 1) && isa(varargin{1},'least_squares_term')
                        lst=varargin{1};
                    else
                        asizec=num2cell(array_size(varargin{1:min(nargin,4)}));
                        lst(asizec{:})=least_squares_term;
                        clear asizec

                        if ischar(varargin{1})
                            varargin{1}=cellstr(varargin{1});
                        end
                        for i=1:numel(lst)
                            if (numel(varargin{1}) > 1)
                                lst(i).descriptor=varargin{1}{i};
                            else
                                lst(i).descriptor=[char(varargin{1}) string_dim(lst,i,'vector')];
                            end
                        end

                        if (nargin >= 2)
                            if ischar(varargin{2})
                                varargin{2}=cellstr(varargin{2});
                            end
                            for i=1:numel(lst)
                                if (numel(varargin{2}) > 1)
                                    lst(i).scale_type=varargin{2}{i};
                                else
                                    lst(i).scale_type=char(varargin{2});
                                end
                            end
                            if (nargin >= 3)
                                for i=1:numel(lst)
                                    if (numel(varargin{3}) > 1)
                                        lst(i).scale     =varargin{3}(i);
                                    else
                                        lst(i).scale     =varargin{3};
                                    end
                                end
                                if (nargin >= 4)
                                    for i=1:numel(lst)
                                        if (numel(varargin{4}) > 1)
                                            lst(i).weight    =varargin{4}(i);
                                        else
                                            lst(i).weight    =varargin{4};
                                        end
                                    end

                                    if (nargin > 4)
                                        warning('least_squares_term:extra_arg',...
                                            'Extra arguments for object of class ''%s''.',...
                                            class(lst));
                                    end
                                end
                            end
                        end
                    end
            end

        end

        function []=disp(lst)

%  display the object

            disp(sprintf('\n'));
            for i=1:numel(lst)
                disp(sprintf('class ''%s'' object ''%s%s'' = \n',...
                    class(lst),inputname(1),string_dim(lst,i)));
                disp(sprintf('    descriptor: ''%s'''  ,lst(i).descriptor));
                disp(sprintf('    scale_type: ''%s'''  ,lst(i).scale_type));
                disp(sprintf('         scale: %g'      ,lst(i).scale));
                disp(sprintf('        weight: %g\n'    ,lst(i).weight));
            end

        end

        function [desc]  =prop_desc(lst,dstr)
            desc=cell(1,numel(lst));
            for i=1:numel(lst)
                if ~isempty(lst(i).descriptor)
                    desc(i)=cellstr(lst(i).descriptor);
                elseif ~isempty(inputname(1))
                    desc(i)=cellstr([inputname(1) string_dim(lst,i,'vector')]);
                elseif exist('dstr','var')
                    desc(i)=cellstr([dstr         string_dim(lst,i,'vector')]);
                else
                    desc(i)=cellstr(['lst'        string_dim(lst,i,'vector')]);
                end
            end
            desc=allempty(desc);
        end
        function [stype] =prop_stype(lst)
            stype=cell(1,numel(lst));
            for i=1:numel(lst)
                stype(i)=cellstr(lst(i).scale_type);
            end
            stype=allequal(stype,'none');
        end
        function [scale] =prop_scale(lst)
            scale=zeros(1,numel(lst));
            for i=1:numel(lst)
                scale(i)=lst(i).scale;
            end
            scale=allequal(scale,1.);
        end
        function [weight]=prop_weight(lst)
            weight=zeros(1,numel(lst));
            for i=1:numel(lst)
                weight(i)=lst(i).weight;
            end
            weight=allequal(weight,1.);
        end
        function [lower] =prop_lower(lst)
            lower=[];
        end
        function [upper] =prop_upper(lst)
            upper=[];
        end
        function [target]=prop_target(lst)
            target=[];
        end
    end

    methods (Static)
        function [rdesc]=dakota_write(fidi,dresp,rdesc)

%  collect only the responses of the appropriate class

            lst=struc_class(dresp,'least_squares_term');

%  write responses

            [rdesc]=rlist_write(fidi,'least_squares_terms','least_squares_term',lst,rdesc);
        end

        function []=dakota_rlev_write(fidi,dresp,params)
        end
    end
end
