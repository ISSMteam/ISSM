%DEBUG class definition
%
%   Usage:
%      debug=debug();

classdef debug
	properties (SetAccess=public) 
		valgrind = false;
		gprof    = false;
		profiling = false;
	end
	methods
		function self = debug(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
				end
			end % }}}
		function self = setdefaultparameters(self) % {{{

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   debug parameters:'));

			fielddisplay(self,'valgrind','use Valgrind to debug (0 or 1)');
			fielddisplay(self,'gprof','use gnu-profiler to find out where the time is spent');
			fielddisplay(self,'profiling','enables profiling (memory, flops, time)');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'object',self,'fieldname','profiling','format','Boolean');
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejsdouble(fid,[modelname '.debug.valgrind'],self.valgrind);
			writejsdouble(fid,[modelname '.debug.gprof'],self.gprof);
			writejsdouble(fid,[modelname '.debug.profiling'],self.profiling);

		end % }}}
	end
end
