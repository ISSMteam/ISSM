%
%  constructor for the linear_equality_constraint class.
%
%  [lec]=linear_equality_constraint(varargin)
%
%  where the required varargin are:
%    matrix        (double row, variable coefficients, NaN)
%    target        (double vector, target values, 0.)
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
classdef linear_equality_constraint
    properties
        matrix    = NaN;
        target    = 0.;
        scale_type='none';
        scale     = 1.;
    end

    methods
        function [lec]=linear_equality_constraint(varargin)

            switch nargin

%  create a default object

                case 0

%  copy the object

                case 1
                    if isa(varargin{1},'linear_equality_constraint')
                        lec=varargin{1};
                    else
                        error('Object ''%s'' is a ''%s'' class object, not ''%s''.',...
                            inputname(1),class(varargin{1}),'linear_equality_constraint');
                    end

%  create the object from the input

                otherwise
                    if     (size(varargin{1},1) == array_numel(varargin{2:min(nargin,4)}) || ...
                            size(varargin{1},1) == 1)
                        asizec=num2cell(array_size(varargin{2:min(nargin,4)}));
                    elseif (array_numel(varargin{2:min(nargin,4)}) == 1)
                        asizec=num2cell([size(varargin{1},1) 1]);
                    else
                        error('Matrix for object of class ''%s'' has inconsistent number of rows.',...
                              class(lec));
                    end
                    lec(asizec{:})=linear_equality_constraint;
                    clear asizec

                    for i=1:numel(lec)
                        if (size(varargin{1},1) > 1)
                            lec(i).matrix    =varargin{1}(i,:);
                        else
                            lec(i).matrix    =varargin{1};
                        end
                    end

                    if (nargin >= 2)
                        for i=1:numel(lec)
                            if (numel(varargin{2}) > 1)
                                lec(i).target    =varargin{2}(i);
                            else
                                lec(i).target    =varargin{2};
                            end
                        end
                        if (nargin >= 3)
                            if ischar(varargin{3})
                                varargin{3}=cellstr(varargin{3});
                            end
                            for i=1:numel(lec)
                                if (numel(varargin{3}) > 1)
                                    lec(i).scale_type=varargin{3}{i};
                                else
                                    lec(i).scale_type=char(varargin{3});
                                end
                            end
                            if (nargin >= 4)
                                for i=1:numel(lec)
                                    if (numel(varargin{4}) > 1)
                                        lec(i).scale     =varargin{4}(i);
                                    else
                                        lec(i).scale     =varargin{4};
                                    end
                                end

                                if (nargin > 4)
                                    warning('linear_equality_constraint:extra_arg',...
                                        'Extra arguments for object of class ''%s''.',...
                                        class(lec));
                                end
                            end
                        end
                    end
            end
        end

        function []=disp(lec)

%  display the object

            disp(sprintf('\n'));
            for i=1:numel(lec)
                disp(sprintf('class ''%s'' object ''%s%s'' = \n',...
                    class(lec),inputname(1),string_dim(lec,i)));
                disp(sprintf('        matrix: %s'      ,string_vec(lec(i).matrix)));
                disp(sprintf('        target: %g'      ,lec(i).target));
                disp(sprintf('    scale_type: ''%s'''  ,lec(i).scale_type));
                disp(sprintf('         scale: %g\n'    ,lec(i).scale));
            end

        end

        function [matrix]=prop_matrix(lec)
            matrix=zeros(numel(lec),0);
            for i=1:numel(lec)
                matrix(i,1:size(lec(i).matrix,2))=lec(i).matrix(1,:);
            end
        end
        function [lower] =prop_lower(lec)
            lower=[];
        end
        function [upper] =prop_upper(lec)
            upper=[];
        end
        function [target]=prop_target(lec)
            target=zeros(size(lec));
            for i=1:numel(lec)
                target(i)=lec(i).target;
            end
            target=allequal(target,0.);
        end
        function [stype] =prop_stype(lec)
            stype=cell(size(lec));
            for i=1:numel(lec)
                stype(i)=cellstr(lec(i).scale_type);
            end
            stype=allequal(stype,'none');
        end
        function [scale] =prop_scale(lec)
            scale=zeros(size(lec));
            for i=1:numel(lec)
                scale(i)=lec(i).scale;
            end
            scale=allequal(scale,1.);
        end
    end

    methods (Static)
        function []=dakota_write(fidi,dvar)

%  collect only the variables of the appropriate class

            lec=struc_class(dvar,'linear_equality_constraint');

%  write constraints

            lclist_write(fidi,'linear_equality_constraints','linear_equality',lec);
        end
    end
end
