%CALVINGTEST class definition
%  For testing calving laws and coefficients
%   Usage:
%      calvingtest=calvingtest();

classdef calvingtest
	properties (SetAccess=public)
		speedfactor = 1;
		independentrate = 0;
	end
	methods
		function self = calvingtest(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					inputstruct=varargin{1};
					list1 = properties('calvingtest');
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
			self.speedfactor = 1;
			self.independentrate = 0;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			%Early return
			if (~strcmp(solution,'TransientSolution') | md.transient.ismovingfront==0), return; end
			md = checkfield(md,'fieldname','calving.speedfactor','>=',0,'NaN',1,'Inf',1, 'singletimeseries', 1);
			md = checkfield(md,'fieldname','calving.independentrate','NaN',1,'Inf',1, 'singletimeseries', 1);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   Calving test parameters:'));
			fielddisplay(self,'speedfactor','calving rate is proportional to the ice velocity (e.g. speedfactor=1 -> calving front should not move)');
			fielddisplay(self,'independentrate','calving rate is independent of the ice velocity.');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			yts=md.constants.yts;
			WriteData(fid,prefix,'name','md.calving.law','data',8,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','speedfactor','format','DoubleMat','mattype',1,'timeserieslength',2,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','independentrate','format','DoubleMat','mattype',1,'timeserieslength',2,'yts',md.constants.yts,'scale',1./yts);
		end % }}}
	end
end
