%CALVINGVONMISESAD class definition
%
%   Usage:
%      calvingvonmisesAD=calvingvonmisesAD();

classdef calvingvonmisesAD
	properties (SetAccess=public) 
		basin_id							  = NaN;
		num_basins           = 0;
		stress_threshold_groundedice = 0.;
		stress_threshold_floatingice = 0.;
		min_thickness = 0.;
	end
	methods
		function self = calvingvonmises(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					inputstruct=varargin{1};
					list1 = properties('calvingvonmisesAD');
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

			self.basin_id             = NaN;
			num_basins           = 0;
			
			%Default sigma max
			self.stress_threshold_groundedice = 1e6;
			self.stress_threshold_floatingice = 150e3;

			%For now we turn this off by setting the threshold to 0
			self.min_thickness = 0.;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			%Early return
			if (~strcmp(solution,'TransientSolution') | md.transient.ismovingfront==0), return; end

			md = checkfield(md,'fieldname','calving.basin_id','Inf',1,'>=',0,'<=',md.frontalforcings.num_basins,'size',[md.mesh.numberofelements 1]);
			md = checkfield(md,'fieldname','frontalforcings.num_basins','numel',1,'NaN',1,'Inf',1,'>',0);
			md = checkfield(md,'fieldname','calving.stress_threshold_groundedice','>',0,'NaN',1,'Inf',1,'size',[md.calving.num_basins 1]);
			md = checkfield(md,'fieldname','calving.stress_threshold_floatingice','>',0,'NaN',1,'Inf',1,'size',[md.calving.num_basins 1]);
			%md = checkfield(md,'fieldname','calving.stress_threshold_groundedice','>',0,'NaN',1,'Inf',1,'size',[1 1]);
			%md = checkfield(md,'fieldname','calving.stress_threshold_floatingice','>',0,'NaN',1,'Inf',1,'size',[1 1]);
			md = checkfield(md,'fieldname','calving.min_thickness','>=',0,'NaN',1,'Inf',1,'numel',1);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   Calving VonMises (AD) parameters:'));
			fielddisplay(self,'basin_id','basin ID for elements');
			fielddisplay(self,'num_basins', 'number of basins');
			fielddisplay(self,'stress_threshold_groundedice','sigma_max applied to grounded ice only [Pa]');
			fielddisplay(self,'stress_threshold_floatingice','sigma_max applied to floating ice only [Pa]');
			fielddisplay(self,'min_thickness','minimum thickness below which no ice is allowed [m]');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			yts=md.constants.yts;
			WriteData(fid,prefix,'name','md.calving.law','data',11,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','num_basins','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','basin_id','data',self.basin_id-1,'name','md.calving.basin_id','format','IntMat','mattype',2); %0-indexed
			%WriteData(fid,prefix,'object',self,'fieldname','stress_threshold_groundedice','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			%WriteData(fid,prefix,'object',self,'fieldname','stress_threshold_floatingice','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			%WriteData(fid,prefix,'object',self,'fieldname','stress_threshold_groundedice','format','Double');
			%WriteData(fid,prefix,'object',self,'fieldname','stress_threshold_floatingice','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','stress_threshold_groundedice','format','DoubleMat');
			WriteData(fid,prefix,'object',self,'fieldname','stress_threshold_floatingice','format','DoubleMat');
			WriteData(fid,prefix,'object',self,'fieldname','min_thickness','format','Double');
		end % }}}
	end
end
