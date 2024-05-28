%GIA class definition for Caron model (Caron et al, Geophysical Journal International, 2017)
%
%	Usage:
%		giacaron=giacaron();
classdef giacaron
	properties (SetAccess=public) 
		
		%Physical constants
		gravitational_constant= NaN;
		surface_radius = NaN;
		core_mantle_boudary_radius= NaN;
		inner_core_boudary_radius= NaN;
		stress_norm= NaN;
		gravity_norm= NaN;
		radius_norm= NaN;
		
		%Numerical parameters
		allow_layer_deletion=NaN;
		verbose_mode=NaN;
		
		%GIA problem setup
		forcing_type=NaN;
		isincompressible=NaN;
		benchmark_mode=NaN;
		calculate_sea_level=NaN;
		calculate_rotational_feedback=NaN;
		subtract_present_day=NaN;
		ntime=NaN;
		nphi=NaN;

		%Earth model
		numlayers   = NaN;
		radius      = NaN;
		lame_mu     = NaN;
		lame_lambda = NaN;
		issolid     = NaN;
		density     = NaN; 
		viscosity   = NaN; 
		isburger    = NaN; 
		transient_viscosity    = NaN; 
		transient_mu		   = NaN; 
		
		%Ice history
		ice_model_identifier=NaN;
		ice_model_ntime=NaN;
		nice_sheets=NaN;
		
	end
	methods
		function self = extrude(self,md) % {{{
		end % }}}
		function self = giacaron(varargin) % {{{
			switch nargin
				case 0
				otherwise
					options=pairoptions(varargin{:});
					body=getfieldvalue(options,'body');
					if strcmpi(body,'earth'), 
						self.numlayers=getfieldvalue(options,'numlayers',5);
						self.calculate_sea_level=true;
						%[self.radius,self.lame_mu, self.lame_lambda, self.issolid, self.density, ...
						%self.viscosity, self.isburger, self.transient_viscosity, self.transient_mu]=...
						%modelinit(self.numlayers);
					elseif strcmpi(body,'europa'), 
						error('giacaron constructor error message: ''europa'' body not implemented yet!');
					else 
						error('giacaron constructor error message: body not implemented yet!');
					end
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{
			 %Physical constants
			gravitational_constant= 6.67259e-11;
			stress_norm= 1e12;
			gravity_norm= 10;
			radius_norm= 1.0;

			%Numerical parameters
			allow_layer_deletion=true;
			verbose_mode=false;

			%GIA problem setup
			forcing_type=11;
			isincompressible=false;
			benchmark_mode=false;
			calculate_sea_level=false;
			calculate_rotational_feedback=true;
			subtract_present_day=true;

			%Earth model
			isburger    = false(self.numlayers,1); 
			
			%Ice history
			nice_sheets=1;

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ~ismember('GiaAnalysis',analyses), return; end
			% Physical constants			
			md = checkfield(md,'fieldname','gia.gravitational_constant','NaN',1,'Inf',1,'numel',1,'>',0);
			md = checkfield(md,'fieldname','gia.surface_radius','NaN',1,'Inf',1,'numel',1,'>',0);
			md = checkfield(md,'fieldname','gia.core_mantle_boundary_radius','Inf',1,'numel',1,'>',0);
			md = checkfield(md,'fieldname','gia.inner_core_boundary_radius','Inf',1,'numel',1,'>',0);
			md = checkfield(md,'fieldname','gia.radius_norm','NaN',1,'Inf',1,'numel',1,'>',0);
			md = checkfield(md,'fieldname','gia.stress_norm','NaN',1,'Inf',1,'numel',1,'>',0);
			md = checkfield(md,'fieldname','gia.gravity_norm','NaN',1,'Inf',1,'numel',1,'>',0);

			%Numerical parameters
			md = checkfield(md,'fieldname','gia.allow_layer_deletion','values',[0 1]);
			md = checkfield(md,'fieldname','gia.verbose_mode','values',[0 1]);

			%GIA problem setup
			md = checkfield(md,'fieldname','gia.forcing_type','NaN',1,'Inf',1,'numel',1,'>',0, '<=', 12);
			md = checkfield(md,'fieldname','gia.isincompressible','values',[0 1]);
			md = checkfield(md,'fieldname','gia.benchmark_mode','values',[0 1]);
			md = checkfield(md,'fieldname','gia.calculate_sea_level','values',[0 1]);
			md = checkfield(md,'fieldname','gia.calculate_rotational_feedback','values',[0 1]);
			md = checkfield(md,'fieldname','gia.subtract_present_day','values',[0 1]);
			md = checkfield(md,'fieldname','gia.ntime','NaN',1,'Inf',1,'numel',1,'>',0);
			md = checkfield(md,'fieldname','gia.ntheta','NaN',1,'Inf',1,'numel',1,'>',0);
			md = checkfield(md,'fieldname','gia.nphi','numel',1,'values', self.ntheta*2);
		
			%Ice history
			md = checkfield(md,'fieldname','gia.ice_model_identifier', 'stringrow', 1, 'size', [1 3]);
			md = checkfield(md,'fieldname','gia.ice_model_ntime', 'NaN', 1, 'Inf', 1, 'numel', 1);
			md = checkfield(md,'fieldname','gia.nice_sheets', 'NaN', 1, 'Inf', 1, '>', 0, 'numel', 1);
			
			%Earth parameters
			md = checkfield(md,'fieldname','gia.numlayers','NaN',1,'Inf',1,'>',0,'numel',1);
			md = checkfield(md,'fieldname','gia.radius','NaN',1,'Inf',1,'size',[md.gia.numlayers 1],'>',0);
			md = checkfield(md,'fieldname','gia.lame_mu','NaN',1,'Inf',1,'size',[md.gia.numlayers 1],'>',0);
			md = checkfield(md,'fieldname','gia.lame_lambda','NaN',1,'Inf',1,'size',[md.gia.numlayers 1],'>',0);
			md = checkfield(md,'fieldname','gia.issolid','NaN',1,'Inf',1,'size',[md.gia.numlayers 1],'>',0,'<',2);
			md = checkfield(md,'fieldname','gia.density','NaN',1,'Inf',1,'size',[md.gia.numlayers 1],'>',0);
			md = checkfield(md,'fieldname','gia.viscosity','NaN',1,'Inf',1,'size',[md.gia.numlayers 1],'>',0);
			md = checkfield(md,'fieldname','gia.isburger','NaN',1,'Inf',1,'size',[md.gia.numlayers 1],'>',0,'<',2);
			md = checkfield(md,'fieldname','gia.transient_viscosity','NaN',1,'Inf',1,'size',[md.gia.numlayers 1],'>',0);
			md = checkfield(md,'fieldname','gia.transient_mu','NaN',1,'Inf',1,'size',[md.gia.numlayers 1],'>',0);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   giacaron parameters:'));
			
			% Physical constants			
			fielddisplay(self,'gravitational_constant', 'Newton''s Gravitational constant, should be equal to approximately 6.67e-11 [m^3.kg^-1.s^-2]');
			fielddisplay(self,'surface_radius','NaN', 'Planet outer surface radius [m]');
			fielddisplay(self,'core_mantle_boundary_radius', 'Planet core mantle boundary radius (optional) [m]'); 
			fielddisplay(self,'inner_core_boundary_radius','Planet inner core boundary radius (optional) [m]'); 
			fielddisplay(self,'radius_norm', 'length normalization constant [m] (default 1.0)');
			fielddisplay(self,'stress_norm', 'stress normalization constant [Pa] (default 1e12)');
			fielddisplay(self,'gravity_norm', 'gravity normalization constant [m.s^-2] (default 10)');

			%Numerical parameters
			fielddisplay(self,'allow_layer_deletion', 'boolean allowing the migration of the starting integration radius while increasing the spherical harmonic degree  (default true)');
			fielddisplay(self,'verbose_mode', 'boolean allowing the program to write more details on terminal (default false)')

			%GIA problem setup
			fielddisplay(self,'forcing_type','integer indicating the nature and depth of the forcing for the Love number calculation: 1: ICB -- Volumic Potential,  2:  ICB -- Pressure,	3:  ICB -- Loading,    4:  ICB -- Tangential traction, 	5:  CMB -- Volumic Potential,  6:  CMB -- Pressure	7:  CMB -- Loading,  8:  CMB -- Tangential traction,	9: SURF -- Volumic Potential, 10: SURF -- Pressure, 11: SURF -- Loading,  12: SURF -- Tangential traction (default 11)'); 
			fielddisplay(self,'isincompressible', 'boolean approximating the mantle rheology to an incompressible body, sets Lame_lambda to 5e14 Pa (default false)');
			fielddisplay(self,'benchmark_mode', 'boolean to enter benchmark mode, writes a lot of outputs from the midst of the calculation (default false)');
			fielddisplay(self,'calculate_sea_level', 'boolean allowing the sea level equation solving (default false)');
			fielddisplay(self,'calculate_rotational_feedback', 'boolean allowing the calculation of rotational feedback (default true)');
			fielddisplay(self,'subtract_present_day', 'boolean, subtracts the present day signal so the calculation is expressed as the difference with the present-day topography, geoid and sea level (default true)');
			fielddisplay(self,'ntime', 'number of time steps')
			fielddisplay(self,'ntheta', 'size of grid in latitude') 
			fielddisplay(self,'nphi', 'size of grid in longitude (should always be ntheta*2)');

			%Ice history
			fielddisplay(self,'ice_model_identifier', 'string identifier for the ice model (3 characters)')
			fielddisplay(self,'ice_model_ntime', 'number of time steps in the original ice model, is used to interpolate the model on the desired time')
			fielddisplay(self,'nice_sheets', 'number of ice regions to be scaled independently')

			%Earth parameters
			fielddisplay(self,'numlayers','number of layers (default 5)');
			fielddisplay(self,'radius','array describing the radius for each interface (numlayers+1) [m]');
			fielddisplay(self,'lame_mu','array describing the shear modulus for each layers (numlayers) [Pa]');
			fielddisplay(self,'lame_lambda','array describing the lame lambda parameter (numlayers) [Pa]');
			fielddisplay(self,'issolid','array describing whether the layer is solid or liquid (default 1) (numlayers)');
			fielddisplay(self,'density','array describing each layer''s density (numlayers) [kg/m^3]');
			fielddisplay(self,'viscosity','array describing each layer''s viscosity (numlayers) [Pa.s]');
			fielddisplay(self,'isburger','array describing whether we adopt a MaxWell (0) or Burgers (1) rheology (default 0)');
			fielddisplay(self,'transient_viscosity','array describing each layer''s transient viscosity, only for Burgers rheologies  (numlayers) [Pa.s]');
			fielddisplay(self,'transient_mu','array describing each layer''s transient shear modulus, only for Burgers rheologies  (numlayers) [Pa]');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'name','md.gia.model','data',2,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','mantle_viscosity','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','lithosphere_thickness','format','DoubleMat','mattype',1,'scale',10^3); %from km to m
			WriteData(fid,prefix,'object',self,'fieldname','cross_section_shape','format','Integer');

			%Physical constants
			WriteData(fid,prefix,'object',self,'fieldname','gravitational_constant','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','surface_radius ','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','core_mantle_boudary_radius','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','inner_core_boudary_radius','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','stress_norm','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','gravity_norm','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','radius_norm','format','Double');

			%Numerical parameters
			WriteData(fid,prefix,'object',self,'fieldname','allow_layer_deletion','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','verbose_mode','format','Boolean');
		
			%GIA problem setup
			WriteData(fid,prefix,'object',self,'fieldname','forcing_type','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','isincompressible','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','benchmark_mode','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','calculate_sea_level','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','calculate_rotational_feedback','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','subtract_present_day','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','ntime','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','ntheta','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','nphi','format','Integer');

			%Earth model
			WriteData(fid,prefix,'object',self,'fieldname','numlayers','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','radius','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','lame_mu','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','lame_lambda ','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','issolid','format','BooleanMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','density','format','DoubleMat','mattype',1); 
			WriteData(fid,prefix,'object',self,'fieldname','viscosity','format','DoubleMat','mattype',1); 
			WriteData(fid,prefix,'object',self,'fieldname','isburger','format','BooleanMat','mattype',1); 
			WriteData(fid,prefix,'object',self,'fieldname','transient_viscosity','format','DoubleMat','mattype',1); 
			WriteData(fid,prefix,'object',self,'fieldname','transient_mu','format','DoubleMat','mattype',1); 
		
			%Inversion
%			dataset_type_dentifier= NaN
			WriteData(fid,prefix,'object',self,'fieldname','ndata_std','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','ndata_minimum_val ','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','ndata_maximum_val ','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','ndata_total','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','n_forward_models','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','ice_coefficients_solving ','format','String');
			WriteData(fid,prefix,'object',self,'fieldname','ninverse_parameters','format','Integer');
		
			%Ice history
			WriteData(fid,prefix,'object',self,'fieldname','ice_model_identifier','format','String');
			WriteData(fid,prefix,'object',self,'fieldname','ice_model_ntime','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','nice_sheets','format','Integer');
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			%Physical constants
			writejsdouble(fid,[modelname '.gia.gravitational_constant'],self.gravitational_constant);
			writejsdouble(fid,[modelname '.gia.surface_radius'],self.surface_radius);
			writejsdouble(fid,[modelname '.gia.core_mantle_boudary_radius'],self.core_mantle_boudary_radius);
			writejsdouble(fid,[modelname '.gia.inner_core_boudary_radius'],self.inner_core_boudary_radius);
			writejsdouble(fid,[modelname '.gia.stress_norm'],self.stress_norm);
			writejsdouble(fid,[modelname '.gia.gravity_norm'],self.gravity_norm);
			writejsdouble(fid,[modelname '.gia.radius_norm'],self.radius_norm);
		
			%Numerical parameters
			writejsdouble(fid,[modelname '.gia.allow_layer_deletion'],self.allow_layer_deletion);
			writejsdouble(fid,[modelname '.gia.verbose_mode'],self.verbose_mode);
		
			%GIA problem setup
			writejsdouble(fid,[modelname '.gia.forcing_type'],self.forcing_type);
			writejsdouble(fid,[modelname '.gia.isincompressible'],self.isincompressible);
			writejsdouble(fid,[modelname '.gia.benchmark_mode'],self.benchmark_mode);
			writejsdouble(fid,[modelname '.gia.calculate_sea_level'],self.calculate_sea_level);
			writejsdouble(fid,[modelname '.gia.calculate_rotational_feedback'],self.calculate_rotational_feedback);
			writejsdouble(fid,[modelname '.gia.subtract_present_day'],self.subtract_present_day);
			writejsdouble(fid,[modelname '.gia.ntime'],self.ntime);
			writejsdouble(fid,[modelname '.gia.ntheta'],self.ntheta);
			writejsdouble(fid,[modelname '.gia.nphi'],self.nphi);

			%Earth model
			writejsdouble(fid,[modelname '.gia.numlayers'],self.numlayers);
			writejsdouble(fid,[modelname '.gia.radius'],self.radius);
			writejsdouble(fid,[modelname '.gia.lame_mu'],self.lame_mu);
			writejsdouble(fid,[modelname '.gia.lame_lambda'],self.lame_lambda);
			writejsdouble(fid,[modelname '.gia.issolid'],self.issolid);
			writejsdouble(fid,[modelname '.gia.density'],self.density); 
			writejsdouble(fid,[modelname '.gia.viscosity'],self.viscosity); 
			writejsdouble(fid,[modelname '.gia.isburger'],self.isburger); 
			writejsdouble(fid,[modelname '.gia.transient_viscosity'],self.transient_viscosity); 
			writejsdouble(fid,[modelname '.gia.transient_mu'],self.transient_mu); 
		
		%Inversion
%		 dataset_type_dentifier= NaN
%			writejsdouble(fid,[modelname '.gia.ndata_std'],self.ndata_std);
%			writejsdouble(fid,[modelname '.gia.ndata_minimum_val'],self.ndata_minimum_val);
%			writejsdouble(fid,[modelname '.gia.ndata_maximum_val'],self.ndata_maximum_val);
%			writejsdouble(fid,[modelname '.gia.ndata_total'],self.ndata_total);
%			writejsdouble(fid,[modelname '.gia.n_forward_models'],self.n_forward_models);
%			writejsdouble(fid,[modelname '.gia.ice_coefficients_solving'],self.ice_coefficients_solving);
%			writejsdouble(fid,[modelname '.gia.ninverse_parameters'],self.ninverse_parameters);
		
			%Ice history
			writejsdouble(fid,[modelname '.gia.ice_model_identifier'],self.ice_model_identifier);
			writejsdouble(fid,[modelname '.gia.ice_model_ntime'],self.ice_model_ntime);
			writejsdouble(fid,[modelname '.gia.nice_sheets'],self.nice_sheets);
			writejsdouble(fid,[modelname '.gia.mantle_viscosity'],self.mantle_viscosity);
			writejsdouble(fid,[modelname '.gia.lithosphere_thickness'],self.lithosphere_thickness);
			writejsdouble(fid,[modelname '.gia.cross_section_shape'],self.cross_section_shape);
			
		end % }}}
	end
end
