%HISTOGRAM BIN UNCERTAIN class definition
%
%	Usage: 
%   hbu=histogram_bin_uncertain('descriptor',descriptor,'pairs_per_variable',pairs_per_variable,...
%           'abscissas',abscissas,'counts',counts,'partition',partition)
%      
%      where hbu is the histogram_bin_uncertain object returned by the constructor, pairs_per_variable the 
%      size of the histogram, 'abscissas' and 'counts' describe the histogram itself. 
%      If the variable is distributed, then a partition vector can be supplied, which can be a 
%      partition vector over elements or vertices. In which case counts, 
%
%   Example:
%      md.qmu.variables.giaid=histogram_bin_uncertain('descriptor','GiammeModelId','pairs_per_variable',3,...
%           'count',[.6 .4 0],''abscissas',[1 2 3]);
%      md.qmu.variables.surfaceloadid=histogram_bin_uncertain(...
%           'descriptor','distributed_SurfaceloadModelId','pairs_per_variable',[2 3 4],...
%           'counts',{[1 0],[.6 .4 0],[.4 .4 .2 0]},...
%           'abscissas',{[1 2],[1 2 3],[1 2 3 4]},...
%           'partition',partition_vector);

classdef histogram_bin_uncertain
	properties
		descriptor='';
		pairs_per_variable={};
		abscissas = {};
		counts = [];
		partition=[];
	end
	methods
		function self=histogram_bin_uncertain(varargin) % {{{
		
			%recover options:
			options=pairoptions(varargin{:});

			%initialize fields:
			self.descriptor=getfieldvalue(options,'descriptor');
			self.pairs_per_variable=getfieldvalue(options,'pairs_per_variable');
			self.abscissas=getfieldvalue(options,'abscissas');
			self.counts=getfieldvalue(options,'counts');

			%if the variable is distributed,  a partition vector should have been 
			%supplied, and that partition vector should have as many partitions 
			%as the pairs_per_variable, abscissas and counts arrays: 
			if self.isdistributed(),
				self.partition=getfieldvalue(options,'partition');
				npart=qmupart2npart(self.partition);
				if npart~=length(self.pairs_per_variable),
					error(sprintf('histogram_bin_uncertain constructor: for the distributed variable  %s the number of pairs_per_variable arrays (%i) should be equal to the number of partitions (%i)',self.descriptor,length(self.pairs_per_variable),npart));
				end
				if npart~=length(self.abscissas),
					error(sprintf('histogram_bin_uncertain constructor: for the distributed variable  %s the number of abscissas arrays (%i) should be equal to the number of partitions (%i)',self.descriptor,length(self.abscissas),npart));
				end
				if npart~=length(self.counts),
					error(sprintf('histogram_bin_uncertain constructor: for the distributed variable  %s the size of counts (%i) should be equal to the number of partitions (%i)',self.descriptor,length(self.counts),npart));
				end
			end
		end % }}}
		function md=checkconsistency(self,md,solution,analyses) % {{{
		end % }}}
			function disp(self) % {{{
			disp(sprintf('   histogram_bin uncertain variable: '));
			fielddisplay(self,'descriptor','name tag');
			fielddisplay(self,'pairs_per_variable','number of bins in histogram');
			fielddisplay(self,'abscissas','abscissas for histogram');
			fielddisplay(self,'counts','probabilities for histogram');
			if ~isempty(self.partition),
				fielddisplay(self,'partition','partition vector defining where sampling will occur');
			end
		end
		%}}}
		%virtual functions needed by qmu processing algorithms:
		%implemented:
        function [desc]=prop_desc(hbu,dstr) % {{{ 
            desc=cell(1,numel(hbu));
            for i=1:numel(hbu)
                if ~isempty(hbu(i).descriptor)
                    desc(i)=cellstr(hbu(i).descriptor);
                elseif ~isempty(inputname(1))
                    desc(i)=cellstr([inputname(1) string_dim(hbu,i,'vector')]);
                elseif exist('dstr','var')
                    desc(i)=cellstr([dstr         string_dim(hbu,i,'vector')]);
                else
                    desc(i)=cellstr(['hbu'        string_dim(hbu,i,'vector')]);
                end
            end
            desc=allempty(desc);
        end  %}}}
        function [mean]=prop_mean(hbu) % {{{
            mean=[];
        end % }}}
        function [stddev]=prop_stddev(hbu) % {{{
            stddev=[];
        end % }}}
        function [lower]=prop_lower(hbu) % {{{
            lower=[];
        end % }}}
        function [upper]=prop_upper(hbu) % {{{
            upper=[];
        end % }}}
		function [pairs_per_variable] =prop_pairs_per_variable(hbu) % {{{
			pairs_per_variable=[];
            for i=1:numel(hbu)
                pairs_per_variable=[pairs_per_variable hbu(i).pairs_per_variable];
            end
        end % }}}
		function [abscissas]=prop_abscissas(hbu) % {{{
			abscissas=[]; 
			for i=1:numel(hbu)
				abscissas=[abscissas hbu(i).abscissas'];
			end
			abscissas=allequal(abscissas,-Inf);
		end % }}}
   		function [counts] =prop_counts(hbu) % {{{
			counts=[]; 
			for i=1:numel(hbu)
				counts=[counts hbu(i).counts'];
			end
			counts=allequal(counts,-Inf);
        end % }}}
        function [initpt]=prop_initpt(hbu) % {{{
            initpt=[];
        end % }}}
        function [initst]=prop_initst(hbu) % {{{ 
            initst=[];
        end % }}}
        function [stype]=prop_stype(hbu) % {{{
            stype={};
        end % }}}
        function [scale]=prop_scale(hbu) % {{{
            scale=[]; 
        end % }}}
		function scaled=isscaled(self) % {{{
			if strncmp(self.descriptor,'scaled_',7),
				scaled=1;
			else
				scaled=0;
			end
		end % }}}
		function distributed=isdistributed(self) % {{{
			if strncmp(self.descriptor,'distributed_',12),
				distributed=1;
			else
				distributed=0;
			end
		end % }}}
	end
    methods (Static)
        function []=dakota_write(fidi,dvar) % {{{
			% collect only the variables of the appropriate class
			hbu=struc_class(dvar,'histogram_bin_uncertain');

			% write variables
            vlist_write(fidi,'histogram_bin_uncertain','hbu',hbu);
        end % }}}
    end
end
