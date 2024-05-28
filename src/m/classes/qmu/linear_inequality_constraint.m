%
%  constructor for the linear_inequality_constraint class.
%
%  [lic]=linear_inequality_constraint(varargin)
%
%  where the required varargin are:
%    matrix        (double row, variable coefficients, NaN)
%    lower         (double vector, lower bounds, -Inf)
%    upper         (double vector, upper bounds, 0.)
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
classdef linear_inequality_constraint
    properties
        matrix    = NaN;
        lower     =-Inf;
        upper     = 0.;
        scale_type='none';
        scale     = 1.;
    end

    methods
        function [lic]=linear_inequality_constraint(varargin)

            switch nargin

%  create a default object

                case 0

%  copy the object

                case 1
                    if isa(varargin{1},'linear_inequality_constraint')
                        lic=varargin{1};
                    else
                        error('Object ''%s'' is a ''%s'' class object, not ''%s''.',...
                            inputname(1),class(varargin{1}),'linear_inequality_constraint');
                    end

%  not enough arguments

                case 2
                    error('Construction of ''%s'' class object requires at least %d inputs.',...
                        'linear_inequality_constraint',3)

%  create the object from the input

                otherwise
                    if     (size(varargin{1},1) == array_numel(varargin{2:min(nargin,5)}) || ...
                            size(varargin{1},1) == 1)
                        asizec=num2cell(array_size(varargin{2:min(nargin,5)}));
                    elseif (array_numel(varargin{2:min(nargin,5)}) == 1)
                        asizec=num2cell([size(varargin{1},1) 1]);
                    else
                        error('Matrix for object of class ''%s'' has inconsistent number of rows.',...
                              class(lic));
                    end
                    lic(asizec{:})=linear_inequality_constraint;
                    clear asizec

                    for i=1:numel(lic)
                        if (size(varargin{1},1) > 1)
                            lic(i).matrix    =varargin{1}(i,:);
                        else
                            lic(i).matrix    =varargin{1};
                        end
                    end

                    if (nargin >= 2)
                        for i=1:numel(lic)
                            if (numel(varargin{2}) > 1)
                                lic(i).lower     =varargin{2}(i);
                            else
                                lic(i).lower     =varargin{2};
                            end
                        end
                        if (nargin >= 3)
                            for i=1:numel(lic)
                                if (numel(varargin{3}) > 1)
                                    lic(i).upper     =varargin{3}(i);
                                else
                                    lic(i).upper     =varargin{3};
                                end
                            end
                            if (nargin >= 4)
                                if ischar(varargin{4})
                                    varargin{4}=cellstr(varargin{4});
                                end
                                for i=1:numel(lic)
                                    if (numel(varargin{4}) > 1)
                                        lic(i).scale_type=varargin{4}{i};
                                    else
                                        lic(i).scale_type=char(varargin{4});
                                    end
                                end
                                if (nargin >= 5)
                                    for i=1:numel(lic)
                                        if (numel(varargin{5}) > 1)
                                            lic(i).scale     =varargin{5}(i);
                                        else
                                            lic(i).scale     =varargin{5};
                                        end
                                    end

                                    if (nargin > 5)
                                        warning('linear_inequality_constraint:extra_arg',...
                                            'Extra arguments for object of class ''%s''.',...
                                            class(lic));
                                    end
                                end
                            end
                        end
                    end
            end
        end

        function []=disp(lic)

%  display the object

            disp(sprintf('\n'));
            for i=1:numel(lic)
                disp(sprintf('class ''%s'' object ''%s%s'' = \n',...
                    class(lic),inputname(1),string_dim(lic,i)));
                disp(sprintf('        matrix: %s'      ,string_vec(lic(i).matrix)));
                disp(sprintf('         lower: %g'      ,lic(i).lower));
                disp(sprintf('         upper: %g'      ,lic(i).upper));
                disp(sprintf('    scale_type: ''%s'''  ,lic(i).scale_type));
                disp(sprintf('         scale: %g\n'    ,lic(i).scale));
            end

        end

        function [matrix]=prop_matrix(lic)
            matrix=zeros(numel(lic),0);
            for i=1:numel(lic)
                matrix(i,1:size(lic(i).matrix,2))=lic(i).matrix(1,:);
            end
        end
        function [lower] =prop_lower(lic)
            lower=zeros(size(lic));
            for i=1:numel(lic)
                lower(i)=lic(i).lower;
            end
            lower=allequal(lower,-Inf);
        end
        function [upper] =prop_upper(lic)
            upper=zeros(size(lic));
            for i=1:numel(lic)
                upper(i)=lic(i).upper;
            end
            upper=allequal(upper,0.);
        end
        function [target]=prop_target(lic)
            target=[];
        end
        function [stype] =prop_stype(lic)
            stype=cell(size(lic));
            for i=1:numel(lic)
                stype(i)=cellstr(lic(i).scale_type);
            end
            stype=allequal(stype,'none');
        end
        function [scale] =prop_scale(lic)
            scale=zeros(size(lic));
            for i=1:numel(lic)
                scale(i)=lic(i).scale;
            end
            scale=allequal(scale,1.);
        end
    end

    methods (Static)
        function []=dakota_write(fidi,dvar)

%  collect only the variables of the appropriate class

            lic=struc_class(dvar,'linear_inequality_constraint');

%  write constraints

            lclist_write(fidi,'linear_inequality_constraints','linear_inequality',lic);
        end
    end
end
