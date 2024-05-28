%
%  definition for the calibration_function class.
%
%  [cf]=calibration_function(varargin)
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
classdef calibration_function
    properties
        descriptor='';
    end

    methods
        function [cf]=calibration_function(varargin)

            switch nargin

%  create a default object

                case 0

%  copy the object or create the object from the input

                otherwise
                    if  (nargin == 1) && isa(varargin{1},'calibration_function')
                        cf=varargin{1};
							else
								asizec=num2cell(array_size(varargin{1:min(nargin,4)}));
								cf(asizec{:})=calibration_function;
								clear asizec

								if ischar(varargin{1})
									varargin{1}=cellstr(varargin{1});
								end
								for i=1:numel(cf)
									if (numel(varargin{1}) > 1)
										cf(i).descriptor=varargin{1}{i};
									else
										cf(i).descriptor=[char(varargin{1}) string_dim(cf,i,'vector')];
									end
								end
                    end
            end

        end

        function []=disp(cf)

%  display the object

            disp(sprintf('\n'));
            for i=1:numel(cf)
                disp(sprintf('class ''%s'' object ''%s%s'' = \n',...
                    class(cf),inputname(1),string_dim(cf,i)));
                disp(sprintf('    descriptor: ''%s'''  ,cf(i).descriptor));
            end

        end

        function [desc]  =prop_desc(cf,dstr)
            desc=cell(1,numel(cf));
            for i=1:numel(cf)
                if ~isempty(cf(i).descriptor)
                    desc(i)=cellstr(cf(i).descriptor);
                elseif ~isempty(inputname(1))
                    desc(i)=cellstr([inputname(1) string_dim(cf,i,'vector')]);
                elseif exist('dstr','var')
                    desc(i)=cellstr([dstr         string_dim(cf,i,'vector')]);
                else
                    desc(i)=cellstr(['cf'         string_dim(cf,i,'vector')]);
                end
            end
            desc=allempty(desc);
			end
			function [stype] =prop_stype(cf)
				stype={};
			end
			function [weight] =prop_weight(cf)
				weight=[];
			end
			function [lower] =prop_lower(cf)
				lower=[];
			end
			function [upper] =prop_upper(cf)
				upper=[];
			end
			function [target] =prop_target(cf)
				target=[];
			end
			function [scale] =prop_scale(cf)
				scale=[];
			end







    end
    methods (Static)
        function [rdesc]=dakota_write(fidi,dresp,rdesc)

%  collect only the responses of the appropriate class

            cf=struc_class(dresp,'calibration_function');

%  write responses

            [rdesc]=rlist_write(fidi,'calibration_terms','calibration_function',cf,rdesc);
        end

        function []=dakota_rlev_write(fidi,dresp,params)
        end
    end
end
