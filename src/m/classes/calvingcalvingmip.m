%CALVINGCALVINGMIP class definition
%   For calvingMIP laws and coefficients
%   Usage:
%      calvingcalvingmip=calvingcalvingmip();

classdef calvingcalvingmip
	properties (SetAccess=public) 
		min_thickness = 0.;
		experiment = 1;
	end
	methods
		function self = calvingcalvingmip(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					inputstruct=varargin{1};
					list1 = properties('calvingcalvingmip');
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
			%For now we turn this off by setting the threshold to 0
			self.min_thickness = 0.;

			self.experiment = 1;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			%Early return
			if (~strcmp(solution,'TransientSolution') | md.transient.ismovingfront==0), return; end

			md = checkfield(md,'fieldname','calving.min_thickness','>=',0,'NaN',1,'Inf',1,'numel',1);
			md = checkfield(md,'fieldname','calving.experiment','values',[0, 1, 2, 3, 4, 5]);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   CalvingMIP parameters:'));
			fielddisplay(self,'experiment','Experiment in CalvingMIP');
			fielddisplay(self,'min_thickness','minimum thickness below which no ice is allowed [m]');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			yts=md.constants.yts;
			WriteData(fid,prefix,'name','md.calving.law','data',12,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','min_thickness','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','experiment','format','Integer');
		end % }}}
	end
end
