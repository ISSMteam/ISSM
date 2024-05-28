%CALVINGDEV2 class definition
%
%   Usage:
%      calvingdev2=calvingdev2();

classdef calvingdev2
	properties (SetAccess=public) 
		stress_threshold_groundedice = 0.;
		stress_threshold_floatingice = 0.;
		height_above_floatation = 0.;
	end
	methods
		function self = calvingdev2(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					inputstruct=varargin{1};
					list1 = properties('calvingdev2');
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

			%Default sigma max
			self.stress_threshold_groundedice = 1.e6;
			self.stress_threshold_floatingice = 150.e3;
			self.height_above_floatation= 0.;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			%Early return
			if (~strcmp(solution,'TransientSolution') | md.transient.ismovingfront==0), return; end

			md = checkfield(md,'fieldname','calving.stress_threshold_groundedice','>',0,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','calving.stress_threshold_floatingice','>',0,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','calving.height_above_floatation','<=',0);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   Calving Dev2 parameters:'));
			fielddisplay(self,'stress_threshold_groundedice','sigma_max applied to grounded ice only [Pa]');
			fielddisplay(self,'stress_threshold_floatingice','sigma_max applied to floating ice only [Pa]');
			fielddisplay(self,'height_above_floatation','height above floatation that the ice tongue must reach before it can calve [m]');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			yts=md.constants.yts;
			WriteData(fid,prefix,'name','md.calving.law','data',7,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','stress_threshold_groundedice','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','stress_threshold_floatingice','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1, 'yts', md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','height_above_floatation', 'format','Double');
		end % }}}
	end
end
