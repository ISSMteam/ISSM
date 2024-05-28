%
%  constructor for the nonlinear_inequality_constraint class.
%
%  [nic]=nonlinear_inequality_constraint(varargin)
%
%  where the required varargin are:
%    descriptor    (char, description, '')
%    lower         (double, lower bound, -Inf)
%    upper         (double, upper bound, 0.)
%  and the optional varargin and defaults are:
%    scale_type    (char, scaling type, 'none')
%    scale         (double, scaling factor, 1.)
%
%  note that zero arguments constructs a default instance; one
%  argument of the class copies the instance; and three or more
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
classdef nonlinear_inequality_constraint
    properties
        descriptor='';
        lower     =-Inf;
        upper     = 0.;
        scale_type='none';
        scale     = 1.;
    end

    methods
        function [nic]=nonlinear_inequality_constraint(varargin)

            switch nargin

 %  create a default object

                case 0

%  copy the object

                case 1
                    if isa(varargin{1},'nonlinear_inequality_constraint')
                        nic=varargin{1};
                    else
                        error('Object ''%s'' is a ''%s'' class object, not ''%s''.',...
                            inputname(1),class(varargin{1}),'nonlinear_inequality_constraint');
                    end

%  not enough arguments

                case 2
                    error('Construction of ''%s'' class object requires at least %d inputs.',...
                        'nonlinear_inequality_constraint',3)

%  create the object from the input

                otherwise
                    asizec=num2cell(array_size(varargin{1:min(nargin,5)}));
                    nic(asizec{:})=nonlinear_inequality_constraint;
                    clear asizec

                    if ischar(varargin{1})
                        varargin{1}=cellstr(varargin{1});
                    end
                    for i=1:numel(nic)
                        if (numel(varargin{1}) > 1)
                            nic(i).descriptor=varargin{1}{i};
                        else
                            nic(i).descriptor=[char(varargin{1}) string_dim(nic,i,'vector')];
                        end
                        if (numel(varargin{2}) > 1)
                            nic(i).lower     =varargin{2}(i);
                        else
                            nic(i).lower     =varargin{2};
                        end
                        if (numel(varargin{3}) > 1)
                            nic(i).upper     =varargin{3}(i);
                        else
                            nic(i).upper     =varargin{3};
                        end
                    end

                    if (nargin >= 4)
                        if ischar(varargin{4})
                            varargin{4}=cellstr(varargin{4});
                        end
                        for i=1:numel(nic)
                            if (numel(varargin{4}) > 1)
                                nic(i).scale_type=varargin{4}{i};
                            else
                                nic(i).scale_type=char(varargin{4});
                            end
                        end
                        if (nargin >= 5)
                            for i=1:numel(nic)
                                if (numel(varargin{5}) > 1)
                                    nic(i).scale     =varargin{5}(i);
                                else
                                    nic(i).scale     =varargin{5};
                                end
                            end

                            if (nargin > 5)
                                warning('objective_function:extra_arg',...
                                    'Extra arguments for object of class ''%s''.',...
                                    class(nic));
                            end
                        end
                    end
            end

        end

        function []=disp(nic)

%  display the object

            disp(sprintf('\n'));
            for i=1:numel(nic)
                disp(sprintf('class ''%s'' object ''%s%s'' = \n',...
                    class(nic),inputname(1),string_dim(nic,i)));
                disp(sprintf('    descriptor: ''%s'''  ,nic(i).descriptor));
                disp(sprintf('         lower: %g'      ,nic(i).lower));
                disp(sprintf('         upper: %g'      ,nic(i).upper));
                disp(sprintf('    scale_type: ''%s'''  ,nic(i).scale_type));
                disp(sprintf('         scale: %g\n'    ,nic(i).scale));
            end

        end

        function [desc]  =prop_desc(nic,dstr)
            desc=cell(1,numel(nic));
            for i=1:numel(nic)
                if ~isempty(nic(i).descriptor)
                    desc(i)=cellstr(nic(i).descriptor);
                elseif ~isempty(inputname(1))
                    desc(i)=cellstr([inputname(1) string_dim(nic,i,'vector')]);
                elseif exist('dstr','var')
                    desc(i)=cellstr([dstr         string_dim(nic,i,'vector')]);
                else
                    desc(i)=cellstr(['nic'        string_dim(nic,i,'vector')]);
                end
            end
            desc=allempty(desc);
        end
        function [stype] =prop_stype(nic)
            stype=cell(size(nic));
            for i=1:numel(nic)
                stype(i)=cellstr(nic(i).scale_type);
            end
            stype=allequal(stype,'none');
        end
        function [scale] =prop_scale(nic)
            scale=zeros(size(nic));
            for i=1:numel(nic)
                scale(i)=nic(i).scale;
            end
            scale=allequal(scale,1.);
        end
        function [weight]=prop_weight(nic)
            weight=[];
        end
        function [lower] =prop_lower(nic)
            lower=zeros(size(nic));
            for i=1:numel(nic)
                lower(i)=nic(i).lower;
            end
            lower=allequal(lower,-Inf);
        end
        function [upper] =prop_upper(nic)
            upper=zeros(size(nic));
            for i=1:numel(nic)
                upper(i)=nic(i).upper;
            end
            upper=allequal(upper,0.);
        end
        function [target]=prop_target(nic)
            target=[];
        end
    end

    methods (Static)
        function [rdesc]=dakota_write(fidi,dresp,rdesc)

%  collect only the responses of the appropriate class

            nic=struc_class(dresp,'nonlinear_inequality_constraint');

%  write responses

            [rdesc]=rlist_write(fidi,'nonlinear_inequality_constraints','nonlinear_inequality',nic,rdesc);
        end

        function []=dakota_rlev_write(fidi,dresp,params)
        end
    end
end
