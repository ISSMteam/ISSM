%RESPONSE_FUNCTION class definition
%
%   Usage:
%      rf=response_function('descriptor',descriptor,'response_levels',respl,...
%           'probability_levels',probl,'reliability_levels',rell,general_reliability_levels',grell,...
%           'partition',partition);
%      where rf is the response function object returned by the constructor. All options except the 
%      descriptor are optional. A partition can be provided for scaled variables. 
% 
%   Example:
%      md.qmu.responses.maxvel=response_function('descriptor','MaxVel','response_levels',[0],...
%                              'probl',[0.0001 0.001 0.01 0.25 0.5 0.75 0.99 0.999 0.9999]');

classdef response_function
    properties
        descriptor='';
        respl     =[];
        probl     =[];
        rell      =[];
        grell     =[];
        partition = [];
    end
    methods
		function self=response_function(varargin) %constructor {{{

			%recover options: 
			options = pairoptions(varargin{:});

			%initialize fields: 
			self.descriptor=getfieldvalue(options,'descriptor');
			self.respl=getfieldvalue(options,'response_levels',[]);
			self.probl=getfieldvalue(options,'probability_levels',...
				[0.0001 0.001 0.01 0.25 0.5 0.75 0.99 0.999 0.9999]);
			self.rell=getfieldvalue(options,'reliability_levels',[]);
			self.grell=getfieldvalue(options,'general_reliability_levels',[]);
			
			%if the response is scaled, a partition vector should have been supplied.
			if self.isscaled(),
				self.partition=getfieldvalue(options,'partition');
				npart=partition_npart(self.partition);
			end
		end %}}}
        function []=disp(rf)% {{{
            %TODO: Convert the following to fielddisplay
            % display the object
            disp(sprintf('\n'));
            for i=1:numel(rf)
                disp(sprintf('class ''%s'' object ''%s%s'' = \n',...
                    class(rf),inputname(1),string_dim(rf,i)));
                disp(sprintf('    descriptor: ''%s'''  ,rf(i).descriptor));
                disp(sprintf('         respl: %s'      ,string_vec(rf(i).respl)));
                disp(sprintf('         probl: %s'      ,string_vec(rf(i).probl)));
                disp(sprintf('          rell: %s'      ,string_vec(rf(i).rell)));
                disp(sprintf('         grell: %s'    ,string_vec(rf(i).grell)));
                if ~isempty(rf.partition),
					fielddisplay(rf,'partition','partition vector defining where response will be computed');
                end
                disp(sprintf('\n'));
            end

        end % }}}
        function [desc]  =prop_desc(rf,dstr)% {{{

            desc=cell(1,numel(rf));
            for i=1:numel(rf)
                if ~isempty(rf(i).descriptor)
                    desc(i)=cellstr(rf(i).descriptor);
                elseif ~isempty(inputname(1))
                    desc(i)=cellstr([inputname(1) string_dim(rf,i,'vector')]);
                elseif exist('dstr','var')
                    desc(i)=cellstr([dstr         string_dim(rf,i,'vector')]);
                else
                    desc(i)=cellstr(['rf'         string_dim(rf,i,'vector')]);
                end
            end
            desc=allempty(desc);
        end % }}}
        function [stype] =prop_stype(rf)% {{{

            stype={};
        end% }}}
        function [scale] =prop_scale(rf)% {{{

            scale=[];
        end% }}}
        function [weight]=prop_weight(rf) % {{{

            weight=[];
        end% }}}
        function [lower] =prop_lower(rf)% {{{

            lower=[];
        end% }}}
        function [upper] =prop_upper(rf)% {{{

            upper=[];
        end % }}}
        function [target]=prop_target(rf) % {{{
            target=[];
        end % }}}
        function [respl,probl,rell,grell]=prop_levels(rf) % {{{
            respl=cell(1,numel(rf));
            probl=cell(1,numel(rf));
            rell =cell(1,numel(rf));
            grell=cell(1,numel(rf));
            for i=1:numel(rf)
                respl(i)={rf(i).respl};
                probl(i)={rf(i).probl};
                rell (i)={rf(i).rell};
                grell(i)={rf(i).grell};
            end
            respl=allempty(respl);
            probl=allempty(probl);
            rell =allempty(rell);
            grell=allempty(grell);
        end % }}}
		%new methods: 
    	function scaled =isscaled(self) % {{{
    		if strncmp(self.descriptor,'scaled_',7),
    			scaled=1;
    		else
    			scaled=0;
    		end
    	end % }}}
    end
    methods (Static)
        function [rdesc]=dakota_write(fidi,dresp,rdesc) % {{{

			%collect only the responses of the appropriate class
            rf=struc_class(dresp,'response_function');
			%write responses
			[rdesc]=rlist_write(fidi,'response_functions','response_function',rf,rdesc);
		end % }}}
		function []=dakota_rlev_write(fidi,dresp,params) % {{{

				%  collect only the responses of the appropriate class
				rf=struc_class(dresp,'response_function');

				%  write response levels
				rlev_write(fidi,rf,params);
        end % }}}
    end
end
