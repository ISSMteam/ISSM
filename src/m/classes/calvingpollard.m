%CALVINGPOLLARD class definition
%
%   Usage:
%      calvingpollard=calvingpollard();

classdef calvingpollard
	properties (SetAccess=public) 
		rc = 0.;
	end
	methods
		function self = calvingpollard(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					inputstruct=varargin{1};
					list1 = properties('calvingpollard');
					list2 = fieldnames(inputstruct);
					for i=1:length(list1)
						fieldname = list1{i};
						if ismember(fieldname,list2),
							self.(fieldname) = inputstruct.(fieldname);
						end
					end
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = extrude(self,md) % {{{
			return;
		end % }}}
		function self = setdefaultparameters(self) % {{{

			%Default from Pollard and DeConto 2015
			self.rc = .75;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			%Early return
			if (~strcmp(solution,'TransientSolution') | md.transient.ismovingfront==0), return; end

			md = checkfield(md,'fieldname','calving.rc','>',0,'<=',1,'numel',1);
		end % }}}
		function disp(self) % {{
			disp(sprintf('   Calving Pollard and DeConto parameters:'));
			fielddisplay(self,'rc','critical depth/thickness ratio (default is 0.75) ');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'name','md.calving.law','data',10,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','rc','format','Double');
		end % }}}
	end
end
