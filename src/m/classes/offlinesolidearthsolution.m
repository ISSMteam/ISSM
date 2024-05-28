%OFFLINESOLIDEARTHSOLUTION class definition
%
%   Usage:
%      addsol=offlinesolidearthsolution();

classdef offlinesolidearthsolution < solidearthsolution
	properties (SetAccess=public) 
	end
	methods
		function self = offlinesolidearthsolution(varargin) % {{{
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

			if ~ismember('SealevelchangeAnalysis',analyses) | (strcmp(solution,'TransientSolution') & md.solidearth.settings.isgrd==1), 
				error('offlinesolidearthsolution.m::checkconsistency: trying to run GRD patterns while supplying an offline solution for those patterns!'); 
			end
			self.checkconsistency@solidearthsolution(md,solution,analyses);

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   external: offlinesolidearth solution:'));
			self.disp@solidearthsolution();
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			disp(sprintf('   external (offlinesolidearthsolution) solution:'));
			WriteData(fid,prefix,'data',2,'name','md.solidearth.external.nature','format','Integer'); %code 2 for offlinesolidearthsolution  class
			self.marshall@solidearthsolution(prefix,md,fid);
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
			self.savemodeljs@solidearthsolution(fid,modelname);
		end % }}}
		function self = extrude(self,md) % {{{
		end % }}}
	end
end
