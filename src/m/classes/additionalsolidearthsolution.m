%ADDITIONALSOLIDEARTHSOLUTION class definition
%
%   Usage:
%      addsol=additionalsolidearthsolution();

classdef additionalsolidearthsolution < solidearthsolution
	properties (SetAccess=public) 
	end
	methods
		function self = additionalsolidearthsolution(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{
			self.setdefaultparameters@solidearthsolution();
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ~ismember('SealevelchangeAnalysis',analyses) | (strcmp(solution,'TransientSolution') & md.solidearth.settings.isgrd==0), 
				error('additionalsolidearthsolution checkconsistency error message: need to run GRD solution if you are supplying a GRD additional pattern solution');
			end
			self.checkconsistency@solidearthsolution(md,solution,analyses);

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   external: additionalsolidearth solution:'));
			self.disp@solidearthsolution();
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'data',1,'name','md.solidearth.external.nature','format','Integer'); %code 1 for additionalsolidearthsolution class
			self.marshall@solidearthsolution(prefix,md,fid);
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
			self.savemodeljs@solidearthsolution(fid,modelname);
		end % }}}
		function self = extrude(self,md) % {{{
		end % }}}
	end
end
