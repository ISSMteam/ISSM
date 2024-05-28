%
%  constructor for the nonlinear_equality_constraint class.
%
%  [nec]=nonlinear_equality_constraint(varargin)
%
%  where the required varargin are:
%    descriptor    (char, description, '')
%    target        (double, target value, 0.)
%  and the optional varargin and defaults are:
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
classdef nonlinear_equality_constraint
    properties
        descriptor='';
        target    = 0.;
        scale_type='none';
        scale     = 1.;
    end

    methods
        function [nec]=nonlinear_equality_constraint(varargin)

            switch nargin

%  create a default object

                case 0

%  copy the object

                case 1
                    if isa(varargin{1},'nonlinear_equality_constraint')
                        nec=varargin{1};
                    else
                        error('Object ''%s'' is a ''%s'' class object, not ''%s''.',...
                            inputname(1),class(varargin{1}),'nonlinear_equality_constraint');
                    end

%  create the object from the input

                otherwise
                    asizec=num2cell(array_size(varargin{1:min(nargin,4)}));
                    nec(asizec{:})=nonlinear_equality_constraint;
                    clear asizec

                    if ischar(varargin{1})
                        varargin{1}=cellstr(varargin{1});
                    end
                    for i=1:numel(nec)
                        if (numel(varargin{1}) > 1)
                            nec(i).descriptor=varargin{1}{i};
                        else
                            nec(i).descriptor=[char(varargin{1}) string_dim(nec,i,'vector')];
                        end
                        if (numel(varargin{2}) > 1)
                            nec(i).target    =varargin{2}(i);
                        else
                            nec(i).target    =varargin{2};
                        end
                    end

                    if (nargin >= 3)
                        if ischar(varargin{3})
                            varargin{3}=cellstr(varargin{3});
                        end
                        for i=1:numel(nec)
                            if (numel(varargin{3}) > 1)
                                nec(i).scale_type=varargin{3}{i};
                            else
                                nec(i).scale_type=char(varargin{3});
                            end
                        end
                        if (nargin >= 4)
                            for i=1:numel(nec)
                                if (numel(varargin{4}) > 1)
                                    nec(i).scale     =varargin{4}(i);
                                else
                                    nec(i).scale     =varargin{4};
                                end
                            end

                            if (nargin > 4)
                                warning('objective_function:extra_arg',...
                                    'Extra arguments for object of class ''%s''.',...
                                    class(nec));
                            end
                        end
                    end
            end

        end

        function []=disp(nec)

%  display the object

            disp(sprintf('\n'));
            for i=1:numel(nec)
                disp(sprintf('class ''%s'' object ''%s%s'' = \n',...
                    class(nec),inputname(1),string_dim(nec,i)));
                disp(sprintf('    descriptor: ''%s'''  ,nec(i).descriptor));
                disp(sprintf('        target: %g'      ,nec(i).target));
                disp(sprintf('    scale_type: ''%s'''  ,nec(i).scale_type));
                disp(sprintf('         scale: %g\n'    ,nec(i).scale));
            end

        end

        function [desc]  =prop_desc(nec,dstr)
            desc=cell(1,numel(nec));
            for i=1:numel(nec)
                if ~isempty(nec(i).descriptor)
                    desc(i)=cellstr(nec(i).descriptor);
                elseif ~isempty(inputname(1))
                    desc(i)=cellstr([inputname(1) string_dim(nec,i,'vector')]);
                elseif exist('dstr','var')
                    desc(i)=cellstr([dstr         string_dim(nec,i,'vector')]);
                else
                    desc(i)=cellstr(['nec'        string_dim(nec,i,'vector')]);
                end
            end
            desc=allempty(desc);
        end
        function [stype] =prop_stype(nec)
            stype=cell(size(nec));
            for i=1:numel(nec)
                stype(i)=cellstr(nec(i).scale_type);
            end
            stype=allequal(stype,'none');
        end
        function [scale] =prop_scale(nec)
            scale=zeros(size(nec));
            for i=1:numel(nec)
                scale(i)=nec(i).scale;
            end
            scale=allequal(scale,1.);
        end
        function [weight]=prop_weight(nec)
            weight=[];
        end
        function [lower] =prop_lower(nec)
            lower=[];
        end
        function [upper] =prop_upper(nec)
            upper=[];
        end
        function [target]=prop_target(nec)
            target=zeros(size(nec));
            for i=1:numel(nec)
                target(i)=nec(i).target;
            end
            target=allequal(target,0.);
        end
    end

    methods (Static)
        function [rdesc]=dakota_write(fidi,dresp,rdesc)

%  collect only the responses of the appropriate class

            nec=struc_class(dresp,'nonlinear_equality_constraint');

%  write responses

            [rdesc]=rlist_write(fidi,'nonlinear_equality_constraints','nonlinear_equality',nec,rdesc);
        end

        function []=dakota_rlev_write(fidi,dresp,params)
        end
    end
end
