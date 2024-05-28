%PRIVATE class definition
%
%   Usage:
%      private=private();

classdef private
	properties (SetAccess=public) 
		isconsistent = true;
		runtimename  = '';
		bamg         = struct();
		solution     = '';
	end
	methods
		function self = private(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   private parameters: do not change'));

			fielddisplay(self,'isconsistent','is model self consistent');
			fielddisplay(self,'runtimename','name of the run launched');
			fielddisplay(self,'bamg','structure with mesh properties constructed if bamg is used to mesh the domain');
			fielddisplay(self,'solution','type of solution launched');

		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejsdouble(fid,[modelname '.priv.isconsistent'],self.isconsistent);
			writejsstring(fid,[modelname '.priv.runtimename'],self.runtimename);
			writejsstring(fid,[modelname '.priv.solution'],self.solution);
			%no bamg 

		end % }}}
	end
end
