%CALVINGHAB class definition
%
%   Usage:
%      calvinghab=calvinghab();

classdef calvinghab
	properties (SetAccess=public) 
		flotation_fraction = 0.;
	end
	methods
		function self = calvinghab(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					inputstruct=varargin{1};
					list1 = properties('calvinghab');
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

			%fraction q = .15 of the flotation thickness at the terminus by default
			self.flotation_fraction = 0.15;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			%Early return
			if (~strcmp(solution,'TransientSolution') | md.transient.ismovingfront==0), return; end

			md = checkfield(md,'fieldname','calving.flotation_fraction','>=',0,'NaN',1,'Inf',1);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   Calving Pi parameters:'));
			fielddisplay(self,'flotation_fraction','fraction of flotation thickness at the terminus');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			yts=md.constants.yts;
			WriteData(fid,prefix,'name','md.calving.law','data',5,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','flotation_fraction','format','DoubleMat','mattype',1);
		end % }}}
	end
end
