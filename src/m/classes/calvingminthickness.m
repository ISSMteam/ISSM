%CALVINGMINTHICKNESS class definition
%
%   Usage:
%      calvingminthickness=calvingminthickness();

classdef calvingminthickness
	properties (SetAccess=public) 
		min_thickness = 0.;
	end
	methods
		function self = calvingminthickness(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					inputstruct=varargin{1};
					list1 = properties('calvingminthickness');
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
		end % }}}
		function self = setdefaultparameters(self) % {{{

			%minimum thickness is 100 m by default
			self.min_thickness=100.;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			%Early return
			if (~strcmp(solution,'TransientSolution') | md.transient.ismovingfront==0), return; end

			md = checkfield(md,'fieldname','calving.min_thickness','>',0,'NaN',1,'Inf',1,'numel',1);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   Calving Minimum thickness:'));
			fielddisplay(self,'min_thickness','minimum thickness below which no ice is allowed [m]');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			yts=md.constants.yts;
			WriteData(fid,prefix,'name','md.calving.law','data',4,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','min_thickness','format','Double');
		end % }}}
	end
end
