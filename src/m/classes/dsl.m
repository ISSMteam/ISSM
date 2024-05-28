%DSL class definition
%
%   Usage:
%      dsl=dsl(); %dynamic sea level class, based on CMIP5 outputs

classdef dsl
	properties (SetAccess=public) 

		global_average_thermosteric_sea_level; %Corresponds to zostoga field in CMIP5 archives. Specified as a temporally variable quantity (in m).
		sea_surface_height_above_geoid; %Corresponds to zos field in CMIP5 archives. Spatial average is 0. Specified as a spatio-temporally variable quantity (in m).
		sea_water_pressure_at_sea_floor; %Corresponds to bpo field in CMIP5 archives. Specified as a spatio-temporally variable quantity (in m equivalent, not in Pa!).

	end
	methods
		function self = dsl(varargin) % {{{
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

			self.global_average_thermosteric_sea_level=NaN;
			self.sea_surface_height_above_geoid=NaN;
			self.sea_water_pressure_at_sea_floor=NaN;

		end % }}}
		function disp(self) % {{{

			disp(sprintf('   dsl parameters:'));
			fielddisplay(self,'global_average_thermosteric_sea_level','Corresponds to zostoga field in CMIP5 archives. Specified as a temporally variable quantity (in m).');
			fielddisplay(self,'sea_surface_height_above_geoid','Corresponds to zos field in CMIP5 archives. Spatial average is 0. Specified as a spatio-temporally variable quantity (in m).');
			fielddisplay(self,'sea_water_pressure_at_sea_floor','Corresponds to bpo field in CMIP5 archives. Specified as a spatio-temporally variable quantity (in m equivalent, not in Pa!).');

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~ismember('SealevelchangeAnalysis',analyses) | (strcmp(solution,'TransientSolution') & md.transient.isslc == 0) | (md.transient.isoceantransport==0),
				return;
			end
			md = checkfield(md,'fieldname','dsl.global_average_thermosteric_sea_level','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','dsl.sea_surface_height_above_geoid','NaN',1,'Inf',1,'timeseries',1);
			md = checkfield(md,'fieldname','dsl.sea_water_pressure_at_sea_floor','NaN',1,'Inf',1,'timeseries',1);
			
			if md.solidearth.settings.compute_bp_grd==1, 
				md = checkfield(md,'fieldname','dsl.sea_water_pressure_at_sea_floor','empty',1);
			end
			
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			yts=md.constants.yts;
			WriteData(fid,prefix,'name','md.dsl.model','data',1,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','global_average_thermosteric_sea_level','format','DoubleMat','mattype',2,'timeseries',1,'timeserieslength',2,'yts',yts); %mattype 2, because we are sending a GMSL value identical everywhere on each element. 
			WriteData(fid,prefix,'object',self,'fieldname','sea_surface_height_above_geoid','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts); %mattype 1 because we specify DSL at vertex locations.
			WriteData(fid,prefix,'object',self,'fieldname','sea_water_pressure_at_sea_floor','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts); %mattype 1 because we specify bottom pressure at vertex locations.

		end % }}}
		function self = extrude(self,md) % {{{
			self.sea_surface_height_above_geoid=project3d(md,'vector',self.sea_surface_height_above_geoid,'type','node','layer',1);
			self.sea_water_pressure_at_sea_floor=project3d(md,'vector',self.sea_water_pressure_at_sea_floor,'type','node','layer',1);
		end % }}}
		function self = initialize(self,md) % {{{

			if isnan(self.global_average_thermosteric_sea_level)
				self.global_average_thermosteric_sea_level=[0;0];
				disp('      no dsl.global_average_thermosteric_sea_level specified: transient values set to zero');
			end
			if isnan(self.sea_surface_height_above_geoid)
				self.sea_surface_height_above_geoid=[zeros(md.mesh.numberofvertices,1);0];
				disp('      no dsl.sea_surface_height_above_geoid specified: transient values set to zero');
			end
			if isnan(self.sea_water_pressure_at_sea_floor)
				self.sea_water_pressure_at_sea_floor=[zeros(md.mesh.numberofvertices,1);0];
				disp('      no dsl.sea_water_pressure_at_sea_floor specified: transient values set to zero');
			end
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{

			writejs1Darray(fid,[modelname '.dsl.global_average_thermosteric_sea_level'],self.global_average_thermosteric_sea_level);
			writejs1Darray(fid,[modelname '.dsl.sea_surface_height_above_geoid'],self.sea_surface_height_above_geoid);
			writejs1Darray(fid,[modelname '.dsl.sea_water_pressure_at_sea_floor'],self.sea_water_pressure_at_sea_floor);

		end % }}}
	end
end
