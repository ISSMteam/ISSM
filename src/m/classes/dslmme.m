%DSLMME class definition
%
%   Usage:
%      dsl=dslmme(); %dynamic sea level class based on a multi-model ensemble of CMIP5 outputs

classdef dslmme
	properties (SetAccess=public) 

		modelid; %index into the multi-model ensemble, determine which field will be used.
		global_average_thermosteric_sea_level; %Corresponds to zostoga field in CMIP5 archives. Specified as a temporally variable quantity (in m) for each ensemble.
		sea_surface_height_above_geoid; %Corresponds to zos field in CMIP5 archives. Spatial average is 0. Specified as a spatio-temporally variable quantity (in m) for each ensemble.
		sea_water_pressure_at_sea_floor; %Corresponds to bpo field in CMIP5 archives. Specified as a spatio-temporally variable quantity (in m equivalent, not in Pa!) for each ensemble.

	end
	methods
		function self = dslmme(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=structtoobj(dsl(),varargin{1});
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

			self.modelid=0;
			self.global_average_thermosteric_sea_level={};
			self.sea_surface_height_above_geoid={};
			self.sea_water_pressure_at_sea_floor={};

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~ismember('SealevelchangeAnalysis',analyses) | (strcmp(solution,'TransientSolution') & md.transient.isslc == 0) | (md.transient.isoceantransport==0),
				return;
			end
			for i=1:length(self.global_average_thermosteric_sea_level),
				md = checkfield(md,'field',self.global_average_thermosteric_sea_level{i},'NaN',1,'Inf',1);
				md = checkfield(md,'field',self.sea_surface_height_above_geoid{i},'NaN',1,'Inf',1,'timeseries',1);
				md = checkfield(md,'field',self.sea_water_pressure_at_sea_floor{i},'NaN',1,'Inf',1,'timeseries',1);
			end
			md = checkfield(md,'field',self.modelid,'NaN',1,'Inf',1,'>=',1,'<=',length(self.global_average_thermosteric_sea_level));

			if md.solidearth.settings.compute_bp_grd==1, 
				md = checkfield(md,'fieldname','dsl.sea_water_pressure_at_sea_floor','empty',1);
			end

		end % }}}
		function disp(self) % {{{

			disp(sprintf('   dsl mme parameters:'));
			fielddisplay(self,'modelid','index into the multi-model ensemble, determine which field will be used.');
			fielddisplay(self,'global_average_thermosteric_sea_level','Corresponds to zostoga field in CMIP5 archives. Specified as a temporally variable quantity (in m) for each ensemble.');
			fielddisplay(self,'sea_surface_height_above_geoid','Corresponds to zos field in CMIP5 archives. Spatial average is 0. Specified as a spatio-temporally variable quantity (in m) for each ensemble.');
			fielddisplay(self,'sea_water_pressure_at_sea_floor','Corresponds to bpo field in CMIP5 archives. Specified as a spatio-temporally variable quantity (in m equivalent, not in Pa!) for each ensemble.');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			WriteData(fid,prefix,'name','md.dsl.model','data',2,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','modelid','format','Double');
			WriteData(fid,prefix,'name','md.dsl.nummodels','data',length(self.global_average_thermosteric_sea_level),'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','global_average_thermosteric_sea_level','format','MatArray','timeseries',1,'timeserieslength',2);
			WriteData(fid,prefix,'object',self,'fieldname','sea_water_pressure_at_sea_floor','format','MatArray','timeserieslength',md.mesh.numberofvertices+1);
			WriteData(fid,prefix,'object',self,'fieldname','sea_surface_height_above_geoid','format','MatArray','timeserieslength',md.mesh.numberofvertices+1);

		end % }}}
		function self = extrude(self,md) % {{{
			for i=1:length(self.global_average_thermosteric_sea_level),
				self.sea_surface_height_above_geoid{i}=project3d(md,'vector',self.sea_surface_height_above_geoid{i},'type','node','layer',1);
				self.sea_water_pressure_at_sea_floor{i}=project3d(md,'vector',self.sea_water_pressure_at_sea_floor{i},'type','node','layer',1);
			end
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
			
			writejsdouble(fid,[modelname '.dsl.modelid'],self.modelid);
			writejscellarray(fid,[modelname '.dsl.global_average_thermosteric_sea_level'],self.global_average_thermosteric_sea_level);
			writejscellarray(fid,[modelname '.dsl.sea_surface_height_above_geoid'],self.sea_surface_height_above_geoid);
			writejs1Darray(fid,[modelname '.dsl.compute_fingerprints'],self.compute_fingerprints);
			writejscellarray(fid,[modelname '.dsl.sea_water_pressure_at_sea_floor'],self.sea_water_pressure_at_sea_floor);

		end % }}}
	end
end
