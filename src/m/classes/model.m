%MODEL class definition
%
%   Usage:
%      md = model(varargin)

classdef model
	properties (SetAccess=public) %Model fields
		% {{{
		%Careful here: no other class should be used as default value this is a bug of matlab
		mesh             = 0;
		mask             = 0;

		geometry         = 0;
		constants        = 0;
		smb              = 0;
		basalforcings    = 0;
		materials        = 0;
		damage           = 0;
		friction         = 0;
		flowequation     = 0;
		timestepping     = 0;
		initialization   = 0;
		rifts            = 0;
		dsl              = 0;
		solidearth       = 0;

		debug            = 0;
		verbose          = 0;
		settings         = 0;
		toolkits         = 0;
		cluster          = 0;

		balancethickness = 0;
		stressbalance    = 0;
		groundingline    = 0;
		hydrology        = 0;
		debris           = 0;
		masstransport    = 0;
		mmemasstransport = 0;
		thermal          = 0;
		steadystate      = 0;
		transient        = 0;
		levelset         = 0;
		calving          = 0;
		frontalforcings  = 0;
		love             = 0;
		esa              = 0;
		sampling         = 0;

		autodiff         = 0;
		inversion        = 0;
		qmu              = 0;
		amr              = 0;
		results          = 0;
		outputdefinition = 0;
		radaroverlay     = 0;
		miscellaneous    = 0;
		private          = 0;
		stochasticforcing= 0;

		%}}}
	end
	methods (Static)
		function md = loadobj(md) % {{{
			% This function is directly called by matlab when a model object is
			% loaded. If the input is a struct it is an old version of model and
			% old fields must be recovered (make sure they are in the deprecated
			% model properties)

			if verLessThan('matlab','7.9'),
				disp('Warning: your matlab version is old and there is a risk that load does not work correctly');
				disp('         if the model is not loaded correctly, rename temporarily loadobj so that matlab does not use it');

				% This is a Matlab bug: all the fields of md have their default value
				% Example of error message:
				% Warning: Error loading an object of class 'model':
				% Undefined function or method 'exist' for input arguments of type 'cell'
				%
				% This has been fixed in MATLAB 7.9 (R2009b) and later versions
			end

			if isstruct(md)
				disp('Recovering model object from a previous version');
				md = structtomodel(model,md);
			end

			%2012 August 4th
			if isa(md.materials,'materials'),
				disp('Recovering old materials');
				if numel(md.materials.rheology_Z)==1 & isnan(md.materials.rheology_Z),
					md.materials=matice(md.materials);
				else
					md.materials=matdamageice(md.materials);
				end
			end
			%2013 April 12
			if numel(md.stressbalance.loadingforce==1)
				md.stressbalance.loadingforce=0*ones(md.mesh.numberofvertices,3);
			end
			%2013 April 17
			if isa(md.hydrology,'hydrology'),
				disp('Recovering old hydrology class');
				md.hydrology=hydrologyshreve(md.materials);
			end
			%2013 October 9
			if ~isa(md.damage,'damage'),
				md.damage=damage();
				md.damage.D=zeros(md.mesh.numberofvertices,1);
				md.damage.spcdamage=NaN*ones(md.mesh.numberofvertices,1);
			end
			%2013 November 18
			if ~isa(md.outputdefinition,'outputdefinition'),
				md.outputdefinition=outputdefinition();
			end
			%2014 March 26th
			if isa(md.mesh,'mesh'),
				disp('Recovering old mesh class');
				if isprop(md.mesh,'dimension'),
					if md.mesh.dimension==2,
						md.mesh=mesh2d(md.mesh);
					else
						md.mesh=mesh3dprisms(md.mesh);
					end
				else
					md.mesh=mesh2dvertical(md.mesh);
				end
			end
			%2014 November 12
			if isa(md.calving,'double'); md.calving=calving(); end
			%2016 October 11
			if isa(md.esa,'double'); md.esa=esa(); end
			%2017 February 10th
			if isa(md.settings,'settings'), %this 'isa' verification: 2018 October 24th
				if md.settings.solver_residue_threshold==0,
					md.settings.solver_residue_threshold = 1e-6;
				end
			end
			%2017 May 4th
			if isa(md.amr,'double'); md.amr=amr(); end
			%2017 Aug 29th
			if isa(md.love,'double'); md.love=love(); end
			%2017 Oct 26th
			if isa(md.calving,'calvingdev')
				disp('Warning: calvingdev is now calvingvonmises');
				md.calving=calvingvonmises(md.calving); 
			end
			%2017 Dec 21st (needs to be here)
			if isempty(md.settings)
				disp('Warning: md.settings had to be reset, make sure to adjust md.settings.output_frequency and other fields');
				md.settings = issmsettings();
			end
			%2018 Dec 1st
			if md.settings.sb_coupling_frequency==0
				md.settings.sb_coupling_frequency=1;
			end
			%2019 Jan..
			if isa(md.frontalforcings,'double');
				if(isprop('meltingrate',md.calving) & ~isnan(md.calving.meltingrate))
		gia			disp('Warning: md.calving.meltingrate is now in md.frontalforcings');
				end
				md.frontalforcings=frontalforcings(md.calving); 
			end
			%2019 Feb 26
			if isa(md.settings.results_on_nodes,'double')
				if md.settings.results_on_nodes == 0
					md.settings.results_on_nodes = {};
				else
					md.settings.results_on_nodes = {'all'};
				end
			end
			%2019 Mar 28, updated 2021 April 23
			if isa(md.smb,'SMBcomponents') | isa(md.smb,'SMBmeltcomponents') | isa(md.smb,'SMBforcing') | isa(md.smb,'SMBgemb') 
				if any(strcmp(fieldnames(md.smb),'isclimatology'))
					if isa(md.smb.isclimatology,'double')
						if prod(size(md.smb.isclimatology)) ~= 1
							md.smb.isclimatology = 0;
						end
						md.timestepping.cycle_forcing=md.smb.isclimatology;
					end
				end
			end
			%2019 Dec 16
			if isa(md.dsl,'double') 
				md.dsl=dsl();
			end
			%2020 April 24
			if isa(md.smb,'SMBgemb')
				if isa(md.smb.isconstrainsurfaceT,'double')
					if prod(size(md.smb.isconstrainsurfaceT)) ~= 1
						md.smb.isconstrainsurfaceT = 0;
					end
				end
			end
			%2021 February 17
			if isa(md.sampling,'double'); md.sampling=sampling(); end
			%VV
			if ~isa(md.stochasticforcing,'stochasticforcing'); md.stochasticforcing=stochasticforcing(); end
			%2022 Oct 28
			if ~isa(md.debris,'debris'); md.debris=debris(); end
			%Mmetransport: Jun 2022:
			if ~isa(md.mmemasstransport,'mmemasstransport'); md.mmemasstransport=mmemasstransport(); end
			% 2026 February 18
			if isa(md.friction, 'frictionjosh') && md.friction.coefficient_max==0; md.friction.coefficient_max=300; end
			% 2026 February 20
			if isa(md.smb, 'SMBpddSicopolis') && md.smb.pdd_fac_ice==0; md.smb.pdd_fac_ice=7.28; end
			if isa(md.smb, 'SMBpddSicopolis') && md.smb.pdd_fac_snow==0; md.smb.pdd_fac_snow=2.73; end
		end% }}}
	end
	methods
		function md = model(varargin) % {{{

			switch nargin
				case 0
					md=setdefaultparameters(md,'earth');
				otherwise
					options=pairoptions(varargin{:});
					planet=getfieldvalue(options,'planet','earth');
					md=setdefaultparameters(md,planet);
				end

		end
		%}}}
		function disp(self) % {{{
			disp(sprintf('%19s: %-23s -- %s','mesh'            ,['[1x1 ' class(self.mesh) ']'],'mesh properties'));
			disp(sprintf('%19s: %-23s -- %s','mask'            ,['[1x1 ' class(self.mask) ']'],'defines grounded and floating elements'));
			disp(sprintf('%19s: %-23s -- %s','geometry'        ,['[1x1 ' class(self.geometry) ']'],'surface elevation, bedrock topography, ice thickness,...'));
			disp(sprintf('%19s: %-23s -- %s','constants'       ,['[1x1 ' class(self.constants) ']'],'physical constants'));
			disp(sprintf('%19s: %-23s -- %s','smb'             ,['[1x1 ' class(self.smb) ']'],'surface mass balance'));
			disp(sprintf('%19s: %-23s -- %s','basalforcings'   ,['[1x1 ' class(self.basalforcings) ']'],'bed forcings'));
			disp(sprintf('%19s: %-23s -- %s','materials'       ,['[1x1 ' class(self.materials) ']'],'material properties'));
			disp(sprintf('%19s: %-23s -- %s','damage'          ,['[1x1 ' class(self.damage) ']'],'parameters for damage evolution solution'));
			disp(sprintf('%19s: %-23s -- %s','friction'        ,['[1x1 ' class(self.friction) ']'],'basal friction/drag properties'));
			disp(sprintf('%19s: %-23s -- %s','flowequation'    ,['[1x1 ' class(self.flowequation) ']'],'flow equations'));
			disp(sprintf('%19s: %-23s -- %s','timestepping'    ,['[1x1 ' class(self.timestepping) ']'],'time stepping for transient models'));
			disp(sprintf('%19s: %-23s -- %s','initialization'  ,['[1x1 ' class(self.initialization) ']'],'initial guess/state'));
			disp(sprintf('%19s: %-23s -- %s','rifts'           ,['[1x1 ' class(self.rifts) ']'],'rifts properties'));
			disp(sprintf('%19s: %-23s -- %s','solidearth'      ,['[1x1 ' class(self.solidearth) ']'],'solidearth inputs and settings'));
			disp(sprintf('%19s: %-23s -- %s','dsl'             ,['[1x1 ' class(self.dsl) ']'],'dynamic sea-level '));
			disp(sprintf('%19s: %-23s -- %s','debug'           ,['[1x1 ' class(self.debug) ']'],'debugging tools (valgrind, gprof)'));
			disp(sprintf('%19s: %-23s -- %s','verbose'         ,['[1x1 ' class(self.verbose) ']'],'verbosity level in solve'));
			disp(sprintf('%19s: %-23s -- %s','settings'        ,['[1x1 ' class(self.settings) ']'],'settings properties'));
			disp(sprintf('%19s: %-23s -- %s','toolkits'        ,['[1x1 ' class(self.toolkits) ']'],'PETSc options for each solution'));
			disp(sprintf('%19s: %-23s -- %s','cluster'         ,['[1x1 ' class(self.cluster) ']'],'cluster parameters (number of CPUs...)'));
			disp(sprintf('%19s: %-23s -- %s','balancethickness',['[1x1 ' class(self.balancethickness) ']'],'parameters for balancethickness solution'));
			disp(sprintf('%19s: %-23s -- %s','stressbalance'   ,['[1x1 ' class(self.stressbalance) ']'],'parameters for stressbalance solution'));
			disp(sprintf('%19s: %-23s -- %s','groundingline'   ,['[1x1 ' class(self.groundingline) ']'],'parameters for groundingline solution'));
			disp(sprintf('%19s: %-23s -- %s','hydrology'       ,['[1x1 ' class(self.hydrology) ']'],'parameters for hydrology solution'));
			disp(sprintf('%19s: %-23s -- %s','debris' 	   ,['[1x1 ' class(self.debris) ']'],'parameters for debris solution'));
			disp(sprintf('%19s: %-23s -- %s','masstransport'   ,['[1x1 ' class(self.masstransport) ']'],'parameters for masstransport solution'));
			disp(sprintf('%19s: %-23s -- %s','mmemasstransport',['[1x1 ' class(self.mmemasstransport) ']'],'parameters for mmemasstransport solution'));
			disp(sprintf('%19s: %-23s -- %s','thermal'         ,['[1x1 ' class(self.thermal) ']'],'parameters for thermal solution'));
			disp(sprintf('%19s: %-23s -- %s','steadystate'     ,['[1x1 ' class(self.steadystate) ']'],'parameters for steadystate solution'));
			disp(sprintf('%19s: %-23s -- %s','transient'       ,['[1x1 ' class(self.transient) ']'],'parameters for transient solution'));
			disp(sprintf('%19s: %-23s -- %s','levelset'        ,['[1x1 ' class(self.levelset) ']'],'parameters for moving boundaries (level-set method)'));
			disp(sprintf('%19s: %-23s -- %s','calving'         ,['[1x1 ' class(self.calving) ']'],'parameters for calving'));
			disp(sprintf('%19s: %-23s -- %s','frontalforcings' ,['[1x1 ' class(self.frontalforcings) ']'],'parameters for frontalforcings'));
			disp(sprintf('%19s: %-23s -- %s','esa'             ,['[1x1 ' class(self.esa) ']'],'parameters for elastic adjustment solution'));
			disp(sprintf('%19s: %-23s -- %s','love'            ,['[1x1 ' class(self.love) ']'],'parameters for love solution'));
			disp(sprintf('%19s: %-23s -- %s','sampling'        ,['[1x1 ' class(self.sampling) ']'],'parameters for stochastic sampler'));
			disp(sprintf('%19s: %-23s -- %s','autodiff'        ,['[1x1 ' class(self.autodiff) ']'],'automatic differentiation parameters'));
			disp(sprintf('%19s: %-23s -- %s','inversion'       ,['[1x1 ' class(self.inversion) ']'],'parameters for inverse methods'));
			disp(sprintf('%19s: %-23s -- %s','qmu'             ,['[1x1 ' class(self.qmu) ']'],'Dakota properties'));
			disp(sprintf('%19s: %-23s -- %s','amr'             ,['[1x1 ' class(self.amr) ']'],'adaptive mesh refinement properties'));
			disp(sprintf('%19s: %-23s -- %s','outputdefinition',['[1x1 ' class(self.outputdefinition) ']'],'output definition'));
			disp(sprintf('%19s: %-23s -- %s','results'         ,['[1x1 ' class(self.results) ']'],'model results'));
			disp(sprintf('%19s: %-23s -- %s','radaroverlay'    ,['[1x1 ' class(self.radaroverlay) ']'],'radar image for plot overlay'));
			disp(sprintf('%19s: %-23s -- %s','miscellaneous'   ,['[1x1 ' class(self.miscellaneous) ']'],'miscellaneous fields'));
			disp(sprintf('%19s: %-23s -- %s','stochasticforcing',['[1x1 ' class(self.stochasticforcing) ']'],'stochasticity applied to model forcings'));
		end % }}}
		function md = setdefaultparameters(md,planet) % {{{

			%initialize subclasses
			md.mesh             = mesh2d();
			md.mask             = mask();
			md.constants        = constants();
			md.geometry         = geometry();
			md.initialization   = initialization();
			md.smb              = SMBforcing();
			md.basalforcings    = basalforcings();
			md.friction         = friction();
			md.rifts            = rifts();
			md.solidearth       = solidearth(planet);
			md.dsl              = dsl();
			md.timestepping     = timestepping();
			md.groundingline    = groundingline();
			md.materials        = matice();
			md.damage           = damage();
			md.flowequation     = flowequation();
			md.debug            = debug();
			md.verbose          = verbose();
			md.settings         = issmsettings();
			md.toolkits         = toolkits();
			md.cluster          = generic();
			md.balancethickness = balancethickness();
			md.stressbalance    = stressbalance();
			md.hydrology        = hydrologyshreve();
			md.debris           = debris();
			md.masstransport    = masstransport();
			md.mmemasstransport = mmemasstransport();
			md.thermal          = thermal();
			md.steadystate      = steadystate();
			md.transient        = transient();
			md.levelset         = levelset();
			md.calving          = calving();
			md.frontalforcings  = frontalforcings();
			md.love             = love();
			md.esa              = esa();
			md.sampling         = sampling();
			md.autodiff         = autodiff();
			md.inversion        = inversion();
			md.qmu              = qmu();
			md.amr              = amr();
			md.radaroverlay     = radaroverlay();
			md.results          = struct();
			md.outputdefinition = outputdefinition();
			md.miscellaneous    = miscellaneous();
			md.private          = private();
			md.stochasticforcing= stochasticforcing();
		end
		%}}}
		function md = checkmessage(md,string) % {{{
			if(nargout~=1) error('wrong usage, model must be an output'); end
			disp(['model not consistent: ' string]);
			md.private.isconsistent=false;
		end
		%}}}
		function md = collapse(md)% {{{
			%COLLAPSE - collapses a 3d mesh into a 2d mesh
			%
			%   This routine collapses a 3d model into a 2d model
			%   and collapses all the fields of the 3d model by
			%   taking their depth-averaged values
			%
			%   Usage:
			%      md=collapse(md)
			%
			%   See also: EXTRUDE, EXTRACT

			%Check that the model is really a 3d model
			if ~strcmp(md.mesh.elementtype(),'Penta'),
				error('collapse error message: only 3d mesh can be collapsed')
			end

			%Start with changing all the fields from the 3d mesh 

			%dealing with the friction law
			%drag is limited to nodes that are on the bedrock.
			if isa(md.friction,'friction'),
				md.friction.coefficient=project2d(md,md.friction.coefficient,1);
				md.friction.p=project2d(md,md.friction.p,1);
				md.friction.q=project2d(md,md.friction.q,1);
			elseif isa(md.friction,'frictioncoulomb'),
				md.friction.coefficient=project2d(md,md.friction.coefficient,1);
				md.friction.coefficientcoulomb=project2d(md,md.friction.coefficientcoulomb,1);
				md.friction.p=project2d(md,md.friction.p,1);
				md.friction.q=project2d(md,md.friction.q,1);
			elseif isa(md.friction,'frictionhydro'),
				md.friction.q=project2d(md,md.friction.q,1);
				md.friction.C=project2d(md,md.friction.C,1);
				md.friction.As=project2d(md,md.friction.As,1);
				md.friction.effective_pressure=project2d(md,md.friction.effective_pressure,1);
			elseif isa(md.friction,'frictionwaterlayer'),
				md.friction.coefficient=project2d(md,md.friction.coefficient,1);
				md.friction.p=project2d(md,md.friction.p,1);
				md.friction.q=project2d(md,md.friction.q,1);
				md.friction.water_layer=project2d(md,md.friction.water_layer,1);
			elseif isa(md.friction,'frictionweertman'),
				md.friction.C=project2d(md,md.friction.C,1);
				md.friction.m=project2d(md,md.friction.m,1);
			elseif isa(md.friction,'frictionweertmantemp'),
				md.friction.C=project2d(md,md.friction.C,1);
				md.friction.m=project2d(md,md.friction.m,1);
			elseif isa(md.friction,'frictionjosh'),
				md.friction.coefficient=project2d(md,md.friction.coefficient,1);
				md.friction.pressure_adjusted_temperature=project2d(md,md.friction.pressure_adjusted_temperature,1);
			else
				disp('friction type not supported');
			end

			%observations
			if ~isnan(md.inversion.vx_obs),
				md.inversion.vx_obs=project2d(md,md.inversion.vx_obs,md.mesh.numberoflayers);
			end
			if ~isnan(md.inversion.vy_obs),
				md.inversion.vy_obs=project2d(md,md.inversion.vy_obs,md.mesh.numberoflayers);
			end
			if ~isnan(md.inversion.vel_obs),
				md.inversion.vel_obs=project2d(md,md.inversion.vel_obs,md.mesh.numberoflayers);
			end
			if ~isnan(md.inversion.thickness_obs),
				md.inversion.thickness_obs=project2d(md,md.inversion.thickness_obs,md.mesh.numberoflayers);
			end
			if ~isnan(md.inversion.cost_functions_coefficients),
				md.inversion.cost_functions_coefficients=project2d(md,md.inversion.cost_functions_coefficients,md.mesh.numberoflayers);
			end
			if numel(md.inversion.min_parameters)>1,
				md.inversion.min_parameters=project2d(md,md.inversion.min_parameters,md.mesh.numberoflayers);
			end
			if numel(md.inversion.max_parameters)>1,
				md.inversion.max_parameters=project2d(md,md.inversion.max_parameters,md.mesh.numberoflayers);
			end
			if isa(md.smb,'SMBforcing') & ~isnan(md.smb.mass_balance),
				md.smb.mass_balance=project2d(md,md.smb.mass_balance,md.mesh.numberoflayers); 
			elseif isa(md.smb,'SMBhenning') & ~isnan(md.smb.smbref),
				md.smb.smbref=project2d(md,md.smb.smbref,md.mesh.numberoflayers);
			end

			%results
			if ~isnan(md.initialization.vx),
				md.initialization.vx=DepthAverage(md,md.initialization.vx);
			end
			if ~isnan(md.initialization.vy),
				md.initialization.vy=DepthAverage(md,md.initialization.vy);
			end
			if ~isnan(md.initialization.vz),
				md.initialization.vz=DepthAverage(md,md.initialization.vz);
			end
			if ~isnan(md.initialization.vel),
				md.initialization.vel=DepthAverage(md,md.initialization.vel);
			end
			if ~isnan(md.initialization.temperature),
				md.initialization.temperature=DepthAverage(md,md.initialization.temperature);
			end
			if ~isnan(md.initialization.pressure),
				md.initialization.pressure=project2d(md,md.initialization.pressure,1);
			end
			if ~isnan(md.initialization.sediment_head),
				md.initialization.sediment_head=project2d(md,md.initialization.sediment_head,1);
			end
			if ~isnan(md.initialization.epl_head),
				md.initialization.epl_head=project2d(md,md.initialization.epl_head,1);
			end
			if ~isnan(md.initialization.epl_thickness),
				md.initialization.epl_thickness=project2d(md,md.initialization.epl_thickness,1);
			end
			if ~isnan(md.initialization.waterfraction),
				md.initialization.waterfraction=project2d(md,md.initialization.waterfraction,1);
			end
			if ~isnan(md.initialization.watercolumn),
				md.initialization.watercolumn=project2d(md,md.initialization.watercolumn,1);
			end
			if ~isnan(md.initialization.debris),
				md.initialization.debris=project2d(md,md.initialization.debris,1);
			end


			%elementstype
			if ~isnan(md.flowequation.element_equation)
				md.flowequation.element_equation=project2d(md,md.flowequation.element_equation,1);
				md.flowequation.vertex_equation=project2d(md,md.flowequation.vertex_equation,1);
				md.flowequation.borderSSA=project2d(md,md.flowequation.borderSSA,1);
				md.flowequation.borderHO=project2d(md,md.flowequation.borderHO,1);
				md.flowequation.borderFS=project2d(md,md.flowequation.borderFS,1);
			end

			%boundary conditions
			md.stressbalance.spcvx=project2d(md,md.stressbalance.spcvx,md.mesh.numberoflayers);
			md.stressbalance.spcvy=project2d(md,md.stressbalance.spcvy,md.mesh.numberoflayers);
			md.stressbalance.spcvz=project2d(md,md.stressbalance.spcvz,md.mesh.numberoflayers);
			md.stressbalance.referential=project2d(md,md.stressbalance.referential,md.mesh.numberoflayers);
			md.stressbalance.loadingforce=project2d(md,md.stressbalance.loadingforce,md.mesh.numberoflayers);
			if numel(md.masstransport.spcthickness)>1,
				md.masstransport.spcthickness=project2d(md,md.masstransport.spcthickness,md.mesh.numberoflayers);
			end
			if numel(md.damage.spcdamage)>1,
				md.damage.spcdamage=project2d(md,md.damage.spcdamage,md.mesh.numberoflayers);
			end
			if numel(md.levelset.spclevelset)>1,
				md.levelset.spclevelset=project2d(md,md.levelset.spclevelset,md.mesh.numberoflayers);
			end
			md.thermal.spctemperature=project2d(md,md.thermal.spctemperature,md.mesh.numberoflayers);

			% Hydrologydc variables
			if isa(md.hydrology,'hydrologydc');
				md.hydrology.spcsediment_head=project2d(md,md.hydrology.spcsediment_head,1);
				md.hydrology.mask_eplactive_node=project2d(md,md.hydrology.mask_eplactive_node,1);
				md.hydrology.sediment_transmitivity=project2d(md,md.hydrology.sediment_transmitivity,1);
				md.hydrology.basal_moulin_input=project2d(md,md.hydrology.basal_moulin_input,1);
				if(md.hydrology.isefficientlayer==1)
					md.hydrology.spcepl_head=project2d(md,md.hydrology.spcepl_head,1);
				end
			end
			
			%materials
			md.materials.rheology_B=DepthAverage(md,md.materials.rheology_B);
			md.materials.rheology_n=project2d(md,md.materials.rheology_n,1);
			if isprop(md.materials,'rheology_E')
				md.materials.rheology_E=project2d(md,md.materials.rheology_E,1);
			end
			
			%damage: 
			if md.damage.isdamage,
				md.damage.D=DepthAverage(md,md.damage.D);
			end

			%special for thermal modeling:
			if ~isnan(md.basalforcings.groundedice_melting_rate),
				md.basalforcings.groundedice_melting_rate=project2d(md,md.basalforcings.groundedice_melting_rate,1); 
			end
			if isprop(md.basalforcings,'floatingice_melting_rate') & ~isnan(md.basalforcings.floatingice_melting_rate),
				md.basalforcings.floatingice_melting_rate=project2d(md,md.basalforcings.floatingice_melting_rate,1); 
			end
			md.basalforcings.geothermalflux=project2d(md,md.basalforcings.geothermalflux,1); %bedrock only gets geothermal flux

			if isprop(md.calving,'coeff') & ~isnan(md.calving.coeff),
				md.calving.coeff=project2d(md,md.calving.coeff,1); 
			end
			if isprop(md.frontalforcings,'meltingrate') & ~isnan(md.frontalforcings.meltingrate),
				md.frontalforcings.meltingrate=project2d(md,md.frontalforcings.meltingrate,1); 
			end

			%update of connectivity matrix
			md.mesh.average_vertex_connectivity=25;

			%Collapse the mesh
			nodes2d=md.mesh.numberofvertices2d;
			elements2d=md.mesh.numberofelements2d;

			%parameters
			md.geometry.surface=project2d(md,md.geometry.surface,1);
			md.geometry.thickness=project2d(md,md.geometry.thickness,1);
			md.geometry.base=project2d(md,md.geometry.base,1);
			if ~isnan(md.geometry.bed),
				md.geometry.bed=project2d(md,md.geometry.bed,1);
			end
			if ~isnan(md.mask.ocean_levelset),
				md.mask.ocean_levelset=project2d(md,md.mask.ocean_levelset,1);
			end
			if ~isnan(md.mask.ice_levelset),
				md.mask.ice_levelset=project2d(md,md.mask.ice_levelset,1);
			end

			%lat long
			if numel(md.mesh.lat)==md.mesh.numberofvertices,
				md.mesh.lat=project2d(md,md.mesh.lat,1);
			end
			if numel(md.mesh.long)==md.mesh.numberofvertices,
				md.mesh.long=project2d(md,md.mesh.long,1);
			end

			%outputdefinitions
			for i=1:length(md.outputdefinition.definitions)
				if isobject(md.outputdefinition.definitions{i})
					%get subfields
					solutionsubfields=fields(md.outputdefinition.definitions{i});
					for j=1:length(solutionsubfields),
						field=md.outputdefinition.definitions{i}.(solutionsubfields{j});
						if length(field)==md.mesh.numberofvertices | length(field)==md.mesh.numberofelements,
							md.outputdefinition.definitions{i}.(solutionsubfields{j})=project2d(md,md.outputdefinition.definitions{i}.(solutionsubfields{j}),1);
						end
					end
				end
			end

			%Initialize 2d mesh
			mesh=mesh2d();
			mesh.x=md.mesh.x2d;
			mesh.y=md.mesh.y2d;
			mesh.numberofvertices=md.mesh.numberofvertices2d;
			mesh.numberofelements=md.mesh.numberofelements2d;
			mesh.elements=md.mesh.elements2d;
			if numel(md.mesh.lat)==md.mesh.numberofvertices,
				mesh.lat=project2d(md,md.mesh.lat,1);
			end
			if numel(md.mesh.long)==md.mesh.numberofvertices,
				mesh.long=project2d(md,md.mesh.long,1);
			end
			mesh.epsg=md.mesh.epsg;
			if numel(md.mesh.scale_factor)==md.mesh.numberofvertices,
				mesh.scale_factor=project2d(md,md.mesh.scale_factor,1);
			end
			if ~isnan(md.mesh.vertexonboundary),
				mesh.vertexonboundary=project2d(md,md.mesh.vertexonboundary,1);
			end
			if ~isnan(md.mesh.elementconnectivity),
				mesh.elementconnectivity=project2d(md,md.mesh.elementconnectivity,1);
			end
			md.mesh=mesh;
			md.mesh.vertexconnectivity=NodeConnectivity(md.mesh.elements,md.mesh.numberofvertices);
			md.mesh.elementconnectivity=ElementConnectivity(md.mesh.elements,md.mesh.vertexconnectivity);
			md.mesh.segments=contourenvelope(md.mesh);

		end % }}}
		function md2 = extract(md,area,varargin) % {{{
			%extract - extract a model according to an Argus contour or flag list
			%
			%   This routine extracts a submodel from a bigger model with respect to a given contour
			%   md must be followed by the corresponding exp file or flags list
			%   It can either be a domain file (argus type, .exp extension), or an array of element flags. 
			%   If user wants every element outside the domain to be 
			%   extract2d, add '~' to the name of the domain file (ex: '~HO.exp');
			%   an empty string '' will be considered as an empty domain
			%   a string 'all' will be considered as the entire domain
			%
			%   Usage:
			%      md2=extract(md,area);
			%
			%   Examples:
			%      md2=extract(md,'Domain.exp');
			%
			%   See also: EXTRUDE, COLLAPSE

			%copy model
			md1=md;

			%recover optoins: 
			options=pairoptions(varargin{:});

			%some checks
			if ((nargin<2) | (nargout~=1)),
				help extract
				error('extract error message: bad usage');
			end

			%get elements that are inside area
			flag_elem=FlagElements(md1,area);
			if ~any(flag_elem),
				error('extracted model is empty');
			end

			%kick out all elements with 3 dirichlets
			if getfieldvalue(options,'spccheck',1)
				spc_elem=find(~flag_elem);
				spc_node=sort(unique(md1.mesh.elements(spc_elem,:)));
				flag=ones(md1.mesh.numberofvertices,1);
				flag(spc_node)=0;
				pos=find(sum(flag(md1.mesh.elements),2)==0);
				flag_elem(pos)=0;
			end

			%extracted elements and nodes lists
			pos_elem=find(flag_elem);
			pos_node=sort(unique(md1.mesh.elements(pos_elem,:)));

			%keep track of some fields
			numberofvertices1=md1.mesh.numberofvertices;
			numberofelements1=md1.mesh.numberofelements;
			numberofvertices2=length(pos_node);
			numberofelements2=length(pos_elem);
			flag_node=zeros(numberofvertices1,1);
			flag_node(pos_node)=1;

			%Create Pelem and Pnode (transform old nodes in new nodes and same thing for the elements)
			Pelem=zeros(numberofelements1,1);
			Pelem(pos_elem)=[1:numberofelements2]';
			Pnode=zeros(numberofvertices1,1);
			Pnode(pos_node)=[1:numberofvertices2]';

			%renumber the elements (some nodes won't exist anymore)
			elements_1=md1.mesh.elements;
			elements_2=elements_1(pos_elem,:);
			elements_2(:,1)=Pnode(elements_2(:,1));
			elements_2(:,2)=Pnode(elements_2(:,2));
			elements_2(:,3)=Pnode(elements_2(:,3));
			if isa(md1.mesh,'mesh3dprisms'),
				elements_2(:,4)=Pnode(elements_2(:,4));
				elements_2(:,5)=Pnode(elements_2(:,5));
				elements_2(:,6)=Pnode(elements_2(:,6));
			end

			%OK, now create the new model!

			%take every field from model
			md2=md1;

			%automatically modify fields

			%loop over model fields
			model_fields=fields(md1);
			for i=1:length(model_fields),
				%get field
				field=md1.(model_fields{i});
				fieldsize=size(field);
				if isobject(field), %recursive call
					object_fields=fields(md1.(model_fields{i}));
					for j=1:length(object_fields),
						%get field
						field=md1.(model_fields{i}).(object_fields{j});
						fieldsize=size(field);
						%size = number of nodes * n
						if fieldsize(1)==numberofvertices1
							md2.(model_fields{i}).(object_fields{j})=field(pos_node,:);
						elseif (fieldsize(1)==numberofvertices1+1)
							md2.(model_fields{i}).(object_fields{j})=[field(pos_node,:); field(end,:)];
						%size = number of elements * n
						elseif fieldsize(1)==numberofelements1
							md2.(model_fields{i}).(object_fields{j})=field(pos_elem,:);
						elseif (fieldsize(1)==numberofelements1+1)
							md2.(model_fields{i}).(object_fields{j})=[field(pos_elem,:); field(end,:)];
						end
					end
				else
					%size = number of nodes * n
					if fieldsize(1)==numberofvertices1
						md2.(model_fields{i})=field(pos_node,:);
					elseif (fieldsize(1)==numberofvertices1+1)
						md2.(model_fields{i})=[field(pos_node,:); field(end,:)];
					%size = number of elements * n
					elseif fieldsize(1)==numberofelements1
						md2.(model_fields{i})=field(pos_elem,:);
					elseif (fieldsize(1)==numberofelements1+1)
						md2.(model_fields{i})=[field(pos_elem,:); field(end,:)];
					end
				end
			end

			%modify some specific fields

			%Mesh
			md2.mesh.numberofelements=numberofelements2;
			md2.mesh.numberofvertices=numberofvertices2;
			md2.mesh.elements=elements_2;

			%mesh.uppervertex mesh.lowervertex
			if isa(md1.mesh,'mesh3dprisms'),
				md2.mesh.uppervertex=md1.mesh.uppervertex(pos_node);
				pos=find(~isnan(md2.mesh.uppervertex));
				md2.mesh.uppervertex(pos)=Pnode(md2.mesh.uppervertex(pos));

				md2.mesh.lowervertex=md1.mesh.lowervertex(pos_node);
				pos=find(~isnan(md2.mesh.lowervertex));
				md2.mesh.lowervertex(pos)=Pnode(md2.mesh.lowervertex(pos));

				md2.mesh.upperelements=md1.mesh.upperelements(pos_elem);
				pos=find(~isnan(md2.mesh.upperelements));
				md2.mesh.upperelements(pos)=Pelem(md2.mesh.upperelements(pos));

				md2.mesh.lowerelements=md1.mesh.lowerelements(pos_elem);
				pos=find(~isnan(md2.mesh.lowerelements));
				md2.mesh.lowerelements(pos)=Pelem(md2.mesh.lowerelements(pos));
			end

			%Initial 2d mesh
			if isa(md1.mesh,'mesh3dprisms'),
				flag_elem_2d=flag_elem(1:md1.mesh.numberofelements2d);
				pos_elem_2d=find(flag_elem_2d);
				flag_node_2d=flag_node(1:md1.mesh.numberofvertices2d);
				pos_node_2d=find(flag_node_2d);

				md2.mesh.numberofelements2d=length(pos_elem_2d);
				md2.mesh.numberofvertices2d=length(pos_node_2d);
				md2.mesh.elements2d=md1.mesh.elements2d(pos_elem_2d,:);
				md2.mesh.elements2d(:,1)=Pnode(md2.mesh.elements2d(:,1));
				md2.mesh.elements2d(:,2)=Pnode(md2.mesh.elements2d(:,2));
				md2.mesh.elements2d(:,3)=Pnode(md2.mesh.elements2d(:,3));

				md2.mesh.x2d=md1.mesh.x(pos_node_2d);
				md2.mesh.y2d=md1.mesh.y(pos_node_2d);
			end

			%Edges
			if(dimension(md.mesh)==2),
				if size(md2.mesh.edges,2)>1, %do not use ~isnan because there are some NaNs...
					%renumber first two columns
					pos=find(md2.mesh.edges(:,4)~=-1);
					md2.mesh.edges(:  ,1)=Pnode(md2.mesh.edges(:,1));
					md2.mesh.edges(:  ,2)=Pnode(md2.mesh.edges(:,2));
					md2.mesh.edges(:  ,3)=Pelem(md2.mesh.edges(:,3));
					md2.mesh.edges(pos,4)=Pelem(md2.mesh.edges(pos,4));
					%remove edges when the 2 vertices are not in the domain.
					md2.mesh.edges=md2.mesh.edges(find(md2.mesh.edges(:,1) & md2.mesh.edges(:,2)),:);
					%Replace all zeros by -1 in the last two columns
					pos=find(md2.mesh.edges(:,3)==0);
					md2.mesh.edges(pos,3)=-1;
					pos=find(md2.mesh.edges(:,4)==0);
					md2.mesh.edges(pos,4)=-1;
					%Invert -1 on the third column with last column (Also invert first two columns!!)
					pos=find(md2.mesh.edges(:,3)==-1);
					md2.mesh.edges(pos,3)=md2.mesh.edges(pos,4);
					md2.mesh.edges(pos,4)=-1;
					values=md2.mesh.edges(pos,2);
					md2.mesh.edges(pos,2)=md2.mesh.edges(pos,1);
					md2.mesh.edges(pos,1)=values;
					%Finally remove edges that do not belong to any element
					pos=find(md2.mesh.edges(:,3)==-1 & md2.mesh.edges(:,4)==-1);
					md2.mesh.edges(pos,:)=[];
				end
			end

			%Penalties
			if ~isnan(md2.stressbalance.vertex_pairing),
				for i=1:size(md1.stressbalance.vertex_pairing,1);
					md2.stressbalance.vertex_pairing(i,:)=Pnode(md1.stressbalance.vertex_pairing(i,:));
				end
				md2.stressbalance.vertex_pairing=md2.stressbalance.vertex_pairing(find(md2.stressbalance.vertex_pairing(:,1)),:);
			end
			if ~isnan(md2.masstransport.vertex_pairing),
				for i=1:size(md1.masstransport.vertex_pairing,1);
					md2.masstransport.vertex_pairing(i,:)=Pnode(md1.masstransport.vertex_pairing(i,:));
				end
				md2.masstransport.vertex_pairing=md2.masstransport.vertex_pairing(find(md2.masstransport.vertex_pairing(:,1)),:);
			end

			%recreate segments
			if isa(md1.mesh,'mesh2d') | isa(md1.mesh','mesh3dsurface'),
				md2.mesh.vertexconnectivity=NodeConnectivity(md2.mesh.elements,md2.mesh.numberofvertices);
				md2.mesh.elementconnectivity=ElementConnectivity(md2.mesh.elements,md2.mesh.vertexconnectivity);
				md2.mesh.segments=contourenvelope(md2.mesh);
				md2.mesh.vertexonboundary=zeros(numberofvertices2,1);
				md2.mesh.vertexonboundary(md2.mesh.segments(:,1:2))=1;
			else
				%First do the connectivity for the contourenvelope in 2d
				md2.mesh.vertexconnectivity=NodeConnectivity(md2.mesh.elements2d,md2.mesh.numberofvertices2d);
				md2.mesh.elementconnectivity=ElementConnectivity(md2.mesh.elements2d,md2.mesh.vertexconnectivity);
				segments=contourenvelope(md2.mesh);
				md2.mesh.vertexonboundary=zeros(numberofvertices2/md2.mesh.numberoflayers,1);
				md2.mesh.vertexonboundary(segments(:,1:2))=1;
				md2.mesh.vertexonboundary=repmat(md2.mesh.vertexonboundary,md2.mesh.numberoflayers,1);
				%Then do it for 3d as usual
				md2.mesh.vertexconnectivity=NodeConnectivity(md2.mesh.elements,md2.mesh.numberofvertices);
				md2.mesh.elementconnectivity=ElementConnectivity(md2.mesh.elements,md2.mesh.vertexconnectivity);
			end

			%Boundary conditions: Dirichlets on new boundary
			%Catch the elements that have not been extracted
			orphans_elem=find(~flag_elem);
			orphans_node=unique(md1.mesh.elements(orphans_elem,:))';
			%Figure out which node are on the boundary between md2 and md1
			nodestoflag1=intersect(orphans_node,pos_node);
			nodestoflag2=Pnode(nodestoflag1);
			if numel(md1.stressbalance.spcvx)>1 & numel(md1.stressbalance.spcvy)>1 & numel(md1.stressbalance.spcvz)>1,
				if isprop(md1.inversion,'vx_obs') & numel(md1.inversion.vx_obs)>1 & numel(md1.inversion.vy_obs)>1
					disp('NOTE: using observed velocities to create constraints along new boundary');
					md2.stressbalance.spcvx(nodestoflag2)=md2.inversion.vx_obs(nodestoflag2); 
					md2.stressbalance.spcvy(nodestoflag2)=md2.inversion.vy_obs(nodestoflag2);
					%MOLHO
					md2.stressbalance.spcvx_base(nodestoflag2)=md2.inversion.vx_obs(nodestoflag2); 
					md2.stressbalance.spcvy_base(nodestoflag2)=md2.inversion.vy_obs(nodestoflag2);
					md2.stressbalance.spcvx_shear(nodestoflag2)=0.;
					md2.stressbalance.spcvy_shear(nodestoflag2)=0.;
				elseif isprop(md1.initialization,'vx') & numel(md1.initialization.vx)>1 & numel(md1.initialization.vy)>1
					disp('NOTE: using initial velocities to create constraints along new boundary');
					md2.stressbalance.spcvx(nodestoflag2)=md2.initialization.vx(nodestoflag2); 
					md2.stressbalance.spcvy(nodestoflag2)=md2.initialization.vy(nodestoflag2);
					%MOLHO
					md2.stressbalance.spcvx_base(nodestoflag2)=md2.initialization.vx(nodestoflag2); 
					md2.stressbalance.spcvy_base(nodestoflag2)=md2.initialization.vy(nodestoflag2);
					md2.stressbalance.spcvx_shear(nodestoflag2)=0.;
					md2.stressbalance.spcvy_shear(nodestoflag2)=0.;
				else
					md2.stressbalance.spcvx(nodestoflag2)=NaN;
					md2.stressbalance.spcvy(nodestoflag2)=NaN;
					warning('Could not set boundary conditions automatically, please set them manually before solve');
				end
				%put 0 for vz
				md2.stressbalance.spcvz(nodestoflag2)=0;
			end
			if ~isnan(md1.thermal.spctemperature),
				md2.thermal.spctemperature(nodestoflag2,1)=1;
			end

			%Results fields
			if isstruct(md1.results),
				md2.results=struct();
				solutionfields=fields(md1.results);
				for i=1:length(solutionfields),
					if isstruct(md1.results.(solutionfields{i}))
						%get subfields
						% loop over time steps
						for p=1:length(md1.results.(solutionfields{i}))
							current = md1.results.(solutionfields{i})(p);
							solutionsubfields=fields(current);
							for j=1:length(solutionsubfields),
							field=md1.results.(solutionfields{i})(p).(solutionsubfields{j});
							if length(field)==numberofvertices1,
								md2.results.(solutionfields{i})(p).(solutionsubfields{j})=field(pos_node);
							elseif length(field)==numberofelements1,
								md2.results.(solutionfields{i})(p).(solutionsubfields{j})=field(pos_elem);
							else
								md2.results.(solutionfields{i})(p).(solutionsubfields{j})=field;
							end
							end
						end
					else
						field=md1.results.(solutionfields{i});
						if length(field)==numberofvertices1,
							md2.results.(solutionfields{i})=field(pos_node);
						elseif length(field)==numberofelements1,
							md2.results.(solutionfields{i})=field(pos_elem);
						else
							md2.results.(solutionfields{i})=field;
						end
					end
				end
			end

			%OutputDefinitions fields
			for i=1:length(md1.outputdefinition.definitions),
				if isobject(md1.outputdefinition.definitions{i})
					%get subfields
					solutionsubfields=fields(md1.outputdefinition.definitions{i});
					for j=1:length(solutionsubfields),
						field=md1.outputdefinition.definitions{i}.(solutionsubfields{j});
						if length(field)==numberofvertices1,
							md2.outputdefinition.definitions{i}.(solutionsubfields{j})=field(pos_node);
						elseif length(field)==numberofelements1,
							md2.outputdefinition.definitions{i}.(solutionsubfields{j})=field(pos_elem);
						elseif size(field,1)==numberofvertices1+1
							md2.outputdefinition.definitions{i}.(solutionsubfields{j})=[field(pos_node,:); field(end,:)];
						end
					end
				end
			end
			
			%independents
			for i=1:length(md1.autodiff.independents)
				independentfield=fields(md1.autodiff.independents{i});
				for j=1:length(independentfield)
					field=md1.autodiff.independents{i}.(independentfield{j});
					if length(field)==numberofvertices1
						md2.autodiff.independents{i}.(independentfield{j})=field(pos_node);
					elseif length(field)==numberofelements1
						md2.autodiff.independents{i}.(independentfield{j})=field(pos_elem);
					end
				end
			end

			%Keep track of pos_node and pos_elem
			md2.mesh.extractedvertices=pos_node;
			md2.mesh.extractedelements=pos_elem;
		end % }}}
		function md2 = refine(md) % {{{
			%refine - split all triangles into 3 to refine the mesh everywhere
			%
			%   This function only works for 2d triangle meshes
			%
			%   Usage:
			%      md2=refine(md);
			%
			%   See also: EXTRUDE, COLLAPSE, EXTRACT

			%Check incoming 
			if ~strcmp(elementtype(md.mesh),'Tria')
				error('not supported for 3d meshes');
			end

			%copy model
			md2=md;

			disp('Getting edges');
			%initialization of some variables
			nbe   = md.mesh.numberofelements;
			nbv   = md.mesh.numberofvertices;
			index = md.mesh.elements;
			elementslist=1:nbe;
			%1: list of edges
			edges=[index(:,[1,2]); index(:,[2,3]); index(:,[3,1])];
			%2: find unique edges
			[edges,I,J]=unique(sort(edges,2),'rows');
			%3: unique edge numbers
			vec=J;
			%4: unique edge numbers in each triangle (2 triangles sharing the same edge will have the same edge number)
			edges_tria=[vec(elementslist+nbe) vec(elementslist+2*nbe) vec(elementslist)];

			% We divide each element as follows
			%
			%                   e2
			%    n1 ------------+------------ n3
			%       \          / \          /
			%        \    1   /   \   3    /
			%         \      /     \      /
			%          \    /   2   \    /
			%           \  /         \  /
			%         e3 +____________\/ e1
			%             \           /
			%              \         /
			%               \   4   /
			%                \     /
			%                 \   /
			%                   n2

			%Create new coordinates
			disp('Remeshing...');
			x_edges = 0.5*(md.mesh.x(edges(:,1)) + md.mesh.x(edges(:,2)));
			y_edges = 0.5*(md.mesh.y(edges(:,1)) + md.mesh.y(edges(:,2)));
			xnew = [md2.mesh.x;x_edges];
			ynew = [md2.mesh.y;y_edges];
			indexnew = [...
				index(:,1)          nbv+edges_tria(:,3) nbv+edges_tria(:,2);...
				nbv+edges_tria(:,2) nbv+edges_tria(:,3) nbv+edges_tria(:,1);...
				nbv+edges_tria(:,2) nbv+edges_tria(:,1) index(:,3);...
				nbv+edges_tria(:,3) index(:,2)          nbv+edges_tria(:,1)];
			%md2.mesh.numberofelements = 4*nbe;
			%md2.mesh.numberofvertices = nbv + size(edges,1);

			%Call Bamg to update other mesh properties
			[bamgmesh_out bamggeom_out]=BamgConvertMesh(indexnew,xnew,ynew);
			md2.mesh.x              = bamgmesh_out.Vertices(:,1);
			md2.mesh.y              = bamgmesh_out.Vertices(:,2);
			md2.mesh.elements       = bamgmesh_out.Triangles(:,1:3);
			md2.mesh.edges          = bamgmesh_out.IssmEdges;
			md2.mesh.segments       = bamgmesh_out.IssmSegments(:,1:3);
			md2.mesh.segmentmarkers = bamgmesh_out.IssmSegments(:,4);
			md2.mesh.numberofelements = size(md2.mesh.elements,1);
			md2.mesh.numberofvertices = length(md2.mesh.x);
			md2.mesh.numberofedges    = size(md2.mesh.edges,1);
			md2.mesh.vertexonboundary = zeros(md2.mesh.numberofvertices,1);
			md2.mesh.vertexonboundary(md2.mesh.segments(:,1:2)) = 1;

			%Deal with boundary
			md2.mesh.vertexonboundary = [md.mesh.vertexonboundary;sum(md.mesh.vertexonboundary(edges),2)==2];
			md2.mesh.elementconnectivity=bamgmesh_out.ElementConnectivity;
			md2.mesh.elementconnectivity(find(isnan(md2.mesh.elementconnectivity)))=0;
			disp(['   Old number of elements: ' num2str(nbe)]);
			disp(['   New number of elements: ' num2str(4*nbe)]);

			disp('Interpolate all fields');
			numberofvertices1 = md.mesh.numberofvertices;
			numberofelements1 = md.mesh.numberofelements;
			nbv2 = md2.mesh.numberofvertices;

			%Create transformation vectors
			nbedges = size(edges,1);
			i = 1:4*nbe;
			j = repmat([1:nbe],1,4);
			v = ones(4*nbe,1);
			m = 4*nbe;
			n = nbe;
			Pelem = sparse(i,j,v,m,n);
			i = [1:nbv,repmat([nbv+1:nbv+nbedges],1,2)];
			j = [1:nbv edges(:)'];
			v = [ones(nbv,1);1/2*ones(2*nbedges,1)];
			m = md2.mesh.numberofvertices;
			n = nbv;
			Pnode = sparse(i,j,v,m,n);

			%Deal with mesh
			if numel(md.mesh.lat)==md.mesh.numberofvertices
				md2.mesh.lat  = Pnode*md.mesh.lat;
				md2.mesh.long = Pnode*md.mesh.long;
			end
			if numel(md.mesh.scale_factor)==md.mesh.numberofvertices
				md2.mesh.scale_factor=Pnode*md.mesh.scale_factor;
			end

			%loop over model fields (except mesh)
			model_fields=setxor(fields(md),{'mesh'});
			for i=1:length(model_fields),
				%get field
				field=md.(model_fields{i});
				fieldsize=size(field);
				if isobject(field), %recursive call
					object_fields=fields(md.(model_fields{i}));
					for j=1:length(object_fields),
						%get field
						field=md.(model_fields{i}).(object_fields{j});
						fieldsize=size(field);
						%size = number of nodes * n
						if fieldsize(1)==numberofvertices1
							md2.(model_fields{i}).(object_fields{j})=Pnode*field;
						elseif (fieldsize(1)==numberofvertices1+1)
							md2.(model_fields{i}).(object_fields{j})=[Pnode*field(1:end-1,:); field(end,:)];
							%size = number of elements * n
						elseif fieldsize(1)==numberofelements1
							md2.(model_fields{i}).(object_fields{j})=Pelem*field;
						elseif (fieldsize(1)==numberofelements1+1)
							md2.(model_fields{i}).(object_fields{j})=[Pelem*field(1:end-1,:); field(end,:)];
						end
					end
				else
					%size = number of nodes * n
					if fieldsize(1)==numberofvertices1
						md2.(model_fields{i})=Pnode*field;
					elseif (fieldsize(1)==numberofvertices1+1)
						md2.(model_fields{i})=[Pnode*field(1:end-1,:); field(end,:)];
					%size = number of elements * n
					elseif fieldsize(1)==numberofelements1
						md2.(model_fields{i})=Pelem*field;
					elseif (fieldsize(1)==numberofelements1+1)
						md2.(model_fields{i})=[Pelem*field(1:end-1,:); field(end,:)];
					end
				end
			end

			%special case: outputdefinitions
			if ~isempty(md.outputdefinition.definitions)
				for i=1:numel(md.outputdefinition.definitions)
					if isa(md.outputdefinition.definitions{i}, 'cfsurfacesquaretransient')
						field = md.outputdefinition.definitions{i}.observations;
						assert(size(field,1)==md.mesh.numberofvertices+1);
						md2.outputdefinition.definitions{i}.observations = [Pnode*field(1:end-1,:); field(end,:)];
						field = md.outputdefinition.definitions{i}.weights;
						assert(size(field,1)==md.mesh.numberofvertices+1);
						md2.outputdefinition.definitions{i}.weights = [Pnode*field(1:end-1,:); field(end,:)];
					elseif isa(md.outputdefinition.definitions{i}, 'cfdragcoeffabsgrad')
						md2.outputdefinition.definitions{i}.weights = Pnode*md.outputdefinition.definitions{i}.weights;
					else
						disp(['skipping md.outputdefinition.definitions{' num2str(i) '} as its class is not yet supported by model.refine']);
						disp('make sure to amend model manually after refinement if this definition is important');
					end
				end
			end

			%special case: independents
			if ~isempty(md.autodiff.independents)
				for i=1:numel(md.autodiff.independents)
					if md.autodiff.independents{i}.nods == md.mesh.numberofvertices
						md2.autodiff.independents{i}.nods = md2.mesh.numberofvertices;
					end
					if numel(md.autodiff.independents{i}.min_parameters)==md.mesh.numberofvertices;
						md2.autodiff.independents{i}.min_parameters = Pnode*md.autodiff.independents{i}.min_parameters;
						md2.autodiff.independents{i}.max_parameters = Pnode*md.autodiff.independents{i}.max_parameters;
					end
				end
			end


		end % }}}
		function md = extrude(md,varargin) % {{{
			%EXTRUDE - vertically extrude a 2d mesh
			%
			%   vertically extrude a 2d mesh and create corresponding 3d mesh.
			%   The vertical distribution can:
			%    - follow a polynomial law
			%    - follow two polynomial laws, one for the lower part and one for the upper part of the mesh
			%    - be discribed by a list of coefficients (between 0 and 1)
			%   
			%
			%   Usage:
			%      md=extrude(md,numlayers,extrusionexponent);
			%      md=extrude(md,numlayers,lowerexponent,upperexponent);
			%      md=extrude(md,listofcoefficients);
			%
			%   Example:
			%      md=extrude(md,15,1.3);
			%      md=extrude(md,15,1.3,1.2);
			%      md=extrude(md,[0 0.2 0.5 0.7 0.9 0.95 1]);
			%
			%   See also: EXTRACT, COLLAPSE

			%some checks on list of arguments
			if ((nargin>4) | (nargin<2) | (nargout~=1)),
				help extrude;
				error('extrude error message');
			end
			if numel(md.geometry.base)~=md.mesh.numberofvertices || numel(md.geometry.surface)~=md.mesh.numberofvertices
				error('model has not been parameterized yet: base and/or surface not set');
			end

			%Extrude the mesh
			if nargin==2, %list of coefficients
				clist=varargin{1};
				if any(clist<0) | any(clist>1),
					error('extrusioncoefficients must be between 0 and 1');
				end
				extrusionlist=sort(unique([clist(:);0;1]));
				numlayers=length(extrusionlist);
			elseif nargin==3, %one polynomial law
				if varargin{2}<=0,
					help extrude;
					error('extrusionexponent must be >=0');
				end
				numlayers=varargin{1};
				extrusionlist=((0:1:numlayers-1)/(numlayers-1)).^varargin{2};
			elseif nargin==4, %two polynomial laws
				numlayers=varargin{1};
				lowerexp=varargin{2};
				upperexp=varargin{3};

				if varargin{2}<=0 | varargin{3}<=0,
					help extrude;
					error('lower and upper extrusionexponents must be >=0');
				end

				lowerextrusionlist=[(0:2/(numlayers-1):1).^lowerexp]/2;
				upperextrusionlist=[(0:2/(numlayers-1):1).^upperexp]/2;
				extrusionlist=sort(unique([lowerextrusionlist 1-upperextrusionlist]));

			end

			if numlayers<2,
				error('number of layers should be at least 2');
			end
			if strcmp(md.mesh.domaintype(),'3D')
				error('Cannot extrude a 3d mesh (extrude cannot be called more than once)');
			end

			%Initialize with 2d mesh
			mesh2d = md.mesh;
			md.mesh=mesh3dprisms();
			md.mesh.x                           = mesh2d.x;
			md.mesh.y                           = mesh2d.y;
			md.mesh.elements                    = mesh2d.elements;
			md.mesh.numberofelements            = mesh2d.numberofelements;
			md.mesh.numberofvertices            = mesh2d.numberofvertices;

			md.mesh.lat                         = mesh2d.lat;
			md.mesh.long                        = mesh2d.long;
			md.mesh.epsg                        = mesh2d.epsg;
			md.mesh.scale_factor                = mesh2d.scale_factor;

			md.mesh.vertexonboundary            = mesh2d.vertexonboundary;
			md.mesh.vertexconnectivity          = mesh2d.vertexconnectivity;
			md.mesh.elementconnectivity         = mesh2d.elementconnectivity;
			md.mesh.average_vertex_connectivity = mesh2d.average_vertex_connectivity;

			md.mesh.extractedvertices           = mesh2d.extractedvertices;
			md.mesh.extractedelements           = mesh2d.extractedelements;

			md.mesh.segments2d                  = mesh2d.segments;

			x3d=[]; 
			y3d=[];
			z3d=[];  %the lower node is on the bed
			thickness3d=md.geometry.thickness; %thickness and bed for these nodes
			bed3d=md.geometry.base;

			%Create the new layers
			for i=1:numlayers,
				x3d=[x3d; md.mesh.x]; 
				y3d=[y3d; md.mesh.y];
				%nodes are distributed between bed and surface accordingly to the given exponent
				z3d=[z3d; bed3d+thickness3d*extrusionlist(i)]; 
			end
			number_nodes3d=size(x3d,1); %number of 3d nodes for the non extruded part of the mesh

			%Extrude elements 
			elements3d=[];
			for i=1:numlayers-1,
				elements3d=[elements3d;[md.mesh.elements+(i-1)*md.mesh.numberofvertices md.mesh.elements+i*md.mesh.numberofvertices]]; %Create the elements of the 3d mesh for the non extruded part
			end
			number_el3d=size(elements3d,1); %number of 3d nodes for the non extruded part of the mesh

			%Keep a trace of lower and upper nodes
			lowervertex=NaN*ones(number_nodes3d,1);
			uppervertex=NaN*ones(number_nodes3d,1);
			lowervertex(md.mesh.numberofvertices+1:end)=1:(numlayers-1)*md.mesh.numberofvertices;
			uppervertex(1:(numlayers-1)*md.mesh.numberofvertices)=md.mesh.numberofvertices+1:number_nodes3d;
			md.mesh.lowervertex=lowervertex;
			md.mesh.uppervertex=uppervertex;

			%same for lower and upper elements
			lowerelements=NaN*ones(number_el3d,1);
			upperelements=NaN*ones(number_el3d,1);
			lowerelements(md.mesh.numberofelements+1:end)=1:(numlayers-2)*md.mesh.numberofelements;
			upperelements(1:(numlayers-2)*md.mesh.numberofelements)=md.mesh.numberofelements+1:(numlayers-1)*md.mesh.numberofelements;
			md.mesh.lowerelements=lowerelements;
			md.mesh.upperelements=upperelements;

			%Save old mesh 
			md.mesh.x2d=md.mesh.x;
			md.mesh.y2d=md.mesh.y;
			md.mesh.elements2d=md.mesh.elements;
			md.mesh.numberofelements2d=md.mesh.numberofelements;
			md.mesh.numberofvertices2d=md.mesh.numberofvertices;

			%Build global 3d mesh 
			md.mesh.elements=elements3d;
			md.mesh.x=x3d;
			md.mesh.y=y3d;
			md.mesh.z=z3d;
			md.mesh.numberofelements=number_el3d;
			md.mesh.numberofvertices=number_nodes3d;
			md.mesh.numberoflayers=numlayers;

			%Ok, now deal with the other fields from the 2d mesh:

			%bedinfo and surface info
			md.mesh.vertexonbase=project3d(md,'vector',ones(md.mesh.numberofvertices2d,1),'type','node','layer',1);
			md.mesh.vertexonsurface=project3d(md,'vector',ones(md.mesh.numberofvertices2d,1),'type','node','layer',md.mesh.numberoflayers);
			md.mesh.vertexonboundary=project3d(md,'vector',md.mesh.vertexonboundary,'type','node');

			%lat long
			md.mesh.lat=project3d(md,'vector',md.mesh.lat,'type','node');
			md.mesh.long=project3d(md,'vector',md.mesh.long,'type','node');
			md.mesh.scale_factor=project3d(md,'vector',md.mesh.scale_factor,'type','node');

			md.geometry=extrude(md.geometry,md);
			md.friction  = extrude(md.friction,md);
			md.inversion = extrude(md.inversion,md);
			md.smb = extrude(md.smb,md);
			md.initialization = extrude(md.initialization,md);

			md.flowequation=md.flowequation.extrude(md);
			md.stressbalance=extrude(md.stressbalance,md);
			md.thermal=md.thermal.extrude(md);
			md.masstransport=md.masstransport.extrude(md);
			md.mmemasstransport=md.mmemasstransport.extrude(md);
			md.levelset=extrude(md.levelset,md);
			md.calving=extrude(md.calving,md);
			md.frontalforcings=extrude(md.frontalforcings,md);
			md.hydrology = extrude(md.hydrology,md);
			md.debris = extrude(md.debris,md);
			md.solidearth = extrude(md.solidearth,md);
			md.dsl = extrude(md.dsl,md);
			md.stochasticforcing = extrude(md.stochasticforcing,md);

			%connectivity
			if ~isnan(md.mesh.elementconnectivity)
				md.mesh.elementconnectivity=repmat(md.mesh.elementconnectivity,numlayers-1,1);
				md.mesh.elementconnectivity(find(md.mesh.elementconnectivity==0))=NaN;
				for i=2:numlayers-1,
					md.mesh.elementconnectivity((i-1)*md.mesh.numberofelements2d+1:(i)*md.mesh.numberofelements2d,:)...
						=md.mesh.elementconnectivity((i-1)*md.mesh.numberofelements2d+1:(i)*md.mesh.numberofelements2d,:)+md.mesh.numberofelements2d;
				end
				md.mesh.elementconnectivity(find(isnan(md.mesh.elementconnectivity)))=0;
			end

			md.materials=extrude(md.materials,md);
			md.damage=extrude(md.damage,md);
			md.mask=extrude(md.mask,md);
			md.qmu=extrude(md.qmu,md);
			md.basalforcings=extrude(md.basalforcings,md);
			md.outputdefinition=extrude(md.outputdefinition,md);

			%increase connectivity if less than 25:
			if md.mesh.average_vertex_connectivity<=25,
				md.mesh.average_vertex_connectivity=100;
			end
		end % }}}
		function md = structtomodel(md,structmd) % {{{

			if ~isstruct(structmd) error('input model is not a structure'); end

			%loaded model is a struct, initialize output and recover all fields
			md = structtoobj(model,structmd);

			%Old field now classes
			if (isfield(structmd,'timestepping') & isnumeric(md.timestepping)), md.timestepping=timestepping(); end
			if (isfield(structmd,'mask') & isnumeric(md.mask)),md.mask=mask(); end

			%Field name change
			if isfield(structmd,'drag'), md.friction.coefficient=structmd.drag; end
			if isfield(structmd,'p'), md.friction.p=structmd.p; end
			if isfield(structmd,'q'), md.friction.q=structmd.p; end
			if isfield(structmd,'melting'), md.basalforcings.floatingice_melting_rate=structmd.melting; end
			if isfield(structmd,'melting_rate'), md.basalforcings.floatingice_melting_rate=structmd.melting_rate; end
			if isfield(structmd,'melting_rate'), md.basalforcings.groundedice_melting_rate=structmd.melting_rate; end
			if isfield(structmd,'accumulation'), md.smb.mass_balance=structmd.accumulation; end
			if isfield(structmd,'numberofgrids'), md.mesh.numberofvertices=structmd.numberofgrids; end
			if isfield(structmd,'numberofgrids2d'), md.mesh.numberofvertices2d=structmd.numberofgrids2d; end
			if isfield(structmd,'uppergrids'), md.mesh.uppervertex=structmd.uppergrids; end
			if isfield(structmd,'lowergrids'), md.mesh.lowervertex=structmd.lowergrids; end
			if isfield(structmd,'gridonbase'), md.mesh.vertexonbase=structmd.gridonbase; end
			if isfield(structmd,'gridonsurface'), md.mesh.vertexonsurface=structmd.gridonsurface; end
			if isfield(structmd,'extractedgrids'), md.mesh.extractedvertices=structmd.extractedgrids; end
			if isfield(structmd,'gridonboundary'), md.mesh.vertexonboundary=structmd.gridonboundary; end
			if isfield(structmd,'petscoptions') & ~isempty(structmd.petscoptions), md.toolkits=structmd.petscoptions; end
			if isfield(structmd,'g'), md.constants.g=structmd.g; end
			if isfield(structmd,'yts'), md.constants.yts=structmd.yts; end
			if isfield(structmd,'surface_mass_balance'), md.smb.mass_balance=structmd.surface_mass_balance; end
			if isfield(structmd,'basal_melting_rate'), md.basalforcings.floatingice_melting_rate=structmd.basal_melting_rate; end
			if isfield(structmd,'geothermalflux'), md.basalforcings.geothermalflux=structmd.geothermalflux; end
			if isfield(structmd,'drag'), md.friction.coefficient=structmd.drag; end
			if isfield(structmd,'drag_coefficient'), md.friction.coefficient=structmd.drag_coefficient; end
			if isfield(structmd,'drag_p'), md.friction.p=structmd.drag_p; end
			if isfield(structmd,'drag_q'), md.friction.q=structmd.drag_q; end
			if isfield(structmd,'riftproperties'), %old implementation
				md.rifts=rifts();
				md.rifts.riftproperties=structmd.riftproperties; 
				md.rifts.riftstruct=structmd.rifts;
				md.rifts.riftproperties=structmd.riftinfo;
			end
			if isfield(structmd,'bamg'), md.private.bamg=structmd.bamg; end
			if isfield(structmd,'lowmem'), md.settings.lowmem=structmd.lowmem; end
			if isfield(structmd,'io_gather'), md.settings.io_gather=structmd.io_gather; end
			if isfield(structmd,'spcwatercolumn'), md.hydrology.spcwatercolumn=structmd.spcwatercolumn; end
			if isfield(structmd,'hydro_n'), md.hydrology.n=structmd.hydro_n; end
			if isfield(structmd,'hydro_p'), md.hydrology.p=structmd.hydro_p; end
			if isfield(structmd,'hydro_q'), md.hydrology.q=structmd.hydro_q; end
			if isfield(structmd,'hydro_CR'), md.hydrology.CR=structmd.hydro_CR; end
			if isfield(structmd,'hydro_kn'), md.hydrology.kn=structmd.hydro_kn; end
			if isfield(structmd,'spctemperature'), md.thermal.spctemperature=structmd.spctemperature; end
			if isfield(structmd,'min_thermal_constraints'), md.thermal.penalty_threshold=structmd.min_thermal_constraints; end
			if isfield(structmd,'artificial_diffusivity'), md.thermal.stabilization=structmd.artificial_diffusivity; end
			if isfield(structmd,'max_nonlinear_iterations'), md.thermal.maxiter=structmd.max_nonlinear_iterations; end
			if isfield(structmd,'stabilize_constraints'), md.thermal.penalty_lock=structmd.stabilize_constraints; end
			if isfield(structmd,'penalty_offset'), md.thermal.penalty_factor=structmd.penalty_offset; end
			if isfield(structmd,'name'), md.miscellaneous.name=structmd.name; end
			if isfield(structmd,'notes'), md.miscellaneous.notes=structmd.notes; end
			if isfield(structmd,'dummy'), md.miscellaneous.dummy=structmd.dummy; end
			if isfield(structmd,'dt'), md.timestepping.time_step=structmd.dt; end
			if isfield(structmd,'ndt'), md.timestepping.final_time=structmd.ndt; end
			if isfield(structmd,'time_adapt'), md.timestepping.time_adapt=structmd.time_adapt; end
			if isfield(structmd,'cfl_coefficient'), md.timestepping.cfl_coefficient=structmd.cfl_coefficient; end
			if isfield(structmd,'spcthickness'), md.masstransport.spcthickness=structmd.spcthickness; end
			if isfield(structmd,'spcthickness'), md.debris.spcthickness=structmd.spcthickness; end
			if isfield(structmd,'artificial_diffusivity'), md.masstransport.stabilization=structmd.artificial_diffusivity; end
			if isfield(structmd,'hydrostatic_adjustment'), md.masstransport.hydrostatic_adjustment=structmd.hydrostatic_adjustment; end
			if isfield(structmd,'penalties'), md.masstransport.vertex_pairing=structmd.penalties; end
			if isfield(structmd,'penalty_offset'), md.masstransport.penalty_factor=structmd.penalty_offset; end
			if isfield(structmd,'deltathickness'), md.mmemasstransport.deltathickness=structmd.deltathickness; end
			if isfield(structmd,'partition'), md.mmemasstransport.partition=structmd.partition; end
			if isfield(structmd,'ids'), md.mmemasstransport.ids=structmd.ids; end
			if isfield(structmd,'B'), md.materials.rheology_B=structmd.B; end
			if isfield(structmd,'n'), md.materials.rheology_n=structmd.n; end
			if isfield(structmd,'rheology_B'), md.materials.rheology_B=structmd.rheology_B; end
			if isfield(structmd,'rheology_n'), md.materials.rheology_n=structmd.rheology_n; end
			if isfield(structmd,'rheology_Z'), md.damage.D=(1-structmd.rheology_Z); end
			if isfield(structmd,'spcthickness'), md.balancethickness.spcthickness=structmd.spcthickness; end
			if isfield(structmd,'artificial_diffusivity'), md.balancethickness.stabilization=structmd.artificial_diffusivity; end
			if isfield(structmd,'dhdt'), md.balancethickness.thickening_rate=structmd.dhdt; end
			if isfield(structmd,'isSIA'), md.flowequation.isSIA=structmd.isSIA; end
			if isfield(structmd,'isFS'), md.flowequation.isFS=structmd.isFS; end
			if isfield(structmd,'elements_type'), md.flowequation.element_equation=structmd.elements_type; end
			if isfield(structmd,'vertices_type'), md.flowequation.vertex_equation=structmd.vertices_type; end
			if isfield(structmd,'eps_rel'), md.steadystate.reltol=structmd.eps_rel; end
			if isfield(structmd,'max_steadystate_iterations'), md.steadystate.maxiter=structmd.max_steadystate_iterations; end
			if isfield(structmd,'isdiagnostic'), md.transient.isstressbalance=structmd.isdiagnostic; end
			if isfield(structmd,'isprognostic'), md.transient.ismasstransport=structmd.isprognostic; end
			if isfield(structmd,'isthermal'), md.transient.isthermal=structmd.isthermal; end
			if isfield(structmd,'control_analysis'), md.inversion.iscontrol=structmd.control_analysis; end
			if isfield(structmd,'weights'), md.inversion.cost_functions_coefficients=structmd.weights; end
			if isfield(structmd,'nsteps'), md.inversion.nsteps=structmd.nsteps; end
			if isfield(structmd,'maxiter_per_step'), md.inversion.maxiter_per_step=structmd.maxiter_per_step; end
			if isfield(structmd,'cm_min'), md.inversion.min_parameters=structmd.cm_min; end
			if isfield(structmd,'cm_max'), md.inversion.max_parameters=structmd.cm_max; end
			if isfield(structmd,'vx_obs'), md.inversion.vx_obs=structmd.vx_obs; end
			if isfield(structmd,'vy_obs'), md.inversion.vy_obs=structmd.vy_obs; end
			if isfield(structmd,'vel_obs'), md.inversion.vel_obs=structmd.vel_obs; end
			if isfield(structmd,'thickness_obs'), md.inversion.thickness_obs=structmd.thickness_obs; end
			if isfield(structmd,'vx'), md.initialization.vx=structmd.vx; end
			if isfield(structmd,'vy'), md.initialization.vy=structmd.vy; end
			if isfield(structmd,'vz'), md.initialization.vz=structmd.vz; end
			if isfield(structmd,'vel'), md.initialization.vel=structmd.vel; end
			if isfield(structmd,'pressure'), md.initialization.pressure=structmd.pressure; end
			if isfield(structmd,'temperature'), md.initialization.temperature=structmd.temperature; end
			if isfield(structmd,'waterfraction'), md.initialization.waterfraction=structmd.waterfraction; end
			if isfield(structmd,'watercolumn'), md.initialization.watercolumn=structmd.watercolumn; end
			if isfield(structmd,'surface'), md.geometry.surface=structmd.surface; end
			if isfield(structmd,'bed'), md.geometry.base=structmd.bed; end
			if isfield(structmd,'thickness'), md.geometry.thickness=structmd.thickness; end
			if isfield(structmd,'bathymetry'), md.geometry.bed=structmd.bathymetry; end
			if isfield(structmd,'thickness_coeff'), md.geometry.hydrostatic_ratio=structmd.thickness_coeff; end
			if isfield(structmd,'connectivity'), md.mesh.average_vertex_connectivity=structmd.connectivity; end
			if isfield(structmd,'extractednodes'), md.mesh.extractedvertices=structmd.extractednodes; end
			if isfield(structmd,'extractedelements'), md.mesh.extractedelements=structmd.extractedelements; end
			if isfield(structmd,'nodeonboundary'), md.mesh.vertexonboundary=structmd.nodeonboundary; end
			if isfield(structmd,'lat'), md.mesh.lat=structmd.lat; end
			if isfield(structmd,'long'), md.mesh.long=structmd.long; end
			if isfield(structmd,'scale_factor'), md.mesh.scale_factor=structmd.scale_factor; end
			if isfield(structmd,'segments'), md.mesh.segments=structmd.segments; end
			if isfield(structmd,'segmentmarkers'), md.mesh.segmentmarkers=structmd.segmentmarkers; end
			if isfield(structmd,'numlayers'), md.mesh.numberoflayers=structmd.numlayers; end
			if isfield(structmd,'numberofelements'), md.mesh.numberofelements=structmd.numberofelements; end
			if isfield(structmd,'numberofvertices'), md.mesh.numberofvertices=structmd.numberofvertices; end
			if isfield(structmd,'numberofnodes'), md.mesh.numberofvertices=structmd.numberofnodes; end
			if isfield(structmd,'numberofedges'), md.mesh.numberofedges=structmd.numberofedges; end
			if isfield(structmd,'numberofelements2d'), md.mesh.numberofelements2d=structmd.numberofelements2d; end
			if isfield(structmd,'numberofnodes2d'), md.mesh.numberofvertices2d=structmd.numberofnodes2d; end
			if isfield(structmd,'nodeconnectivity'), md.mesh.vertexconnectivity=structmd.nodeconnectivity; end
			if isfield(structmd,'elementconnectivity'), md.mesh.elementconnectivity=structmd.elementconnectivity; end
			if isfield(structmd,'uppernodes'), md.mesh.uppervertex=structmd.uppernodes; end
			if isfield(structmd,'lowernodes'), md.mesh.lowervertex=structmd.lowernodes; end
			if isfield(structmd,'upperelements'), md.mesh.upperelements=structmd.upperelements; end
			if isfield(structmd,'lowerelements'), md.mesh.lowerelements=structmd.lowerelements; end
			if isfield(structmd,'nodeonsurface'), md.mesh.vertexonsurface=structmd.nodeonsurface; end
			if isfield(structmd,'nodeonbase'), md.mesh.vertexonbase=structmd.nodeonbase; end
			if isfield(structmd,'elements2d'), md.mesh.elements2d=structmd.elements2d; end
			if isfield(structmd,'y2d'), md.mesh.y2d=structmd.y2d; end
			if isfield(structmd,'x2d'), md.mesh.x2d=structmd.x2d; end
			if isfield(structmd,'elements'), md.mesh.elements=structmd.elements; end
			if isfield(structmd,'edges'), 
				md.mesh.edges=structmd.edges; 
				md.mesh.edges(isnan(md.mesh.edges))=-1;
			end
			if isfield(structmd,'y'), md.mesh.y=structmd.y; end
			if isfield(structmd,'x'), md.mesh.x=structmd.x; end
			if isfield(structmd,'z'), md.mesh.z=structmd.z; end
			if isfield(structmd,'diagnostic_ref'), md.stressbalance.referential=structmd.diagnostic_ref; end
			if isfield(structmd,'npart'); md.qmu.numberofpartitions=structmd.npart; end
			if isfield(structmd,'part'); md.qmu.partition=structmd.part; end

			if isnumeric(md.verbose),
				md.verbose=verbose;
			end

			if isfield(structmd,'spcvelocity'), 
				md.stressbalance.spcvx=NaN*ones(md.mesh.numberofvertices,1);
				md.stressbalance.spcvy=NaN*ones(md.mesh.numberofvertices,1);
				md.stressbalance.spcvz=NaN*ones(md.mesh.numberofvertices,1);
				pos=find(structmd.spcvelocity(:,1)); md.stressbalance.spcvx(pos)=structmd.spcvelocity(pos,4); 
				pos=find(structmd.spcvelocity(:,2)); md.stressbalance.spcvy(pos)=structmd.spcvelocity(pos,5); 
				pos=find(structmd.spcvelocity(:,3)); md.stressbalance.spcvz(pos)=structmd.spcvelocity(pos,6); 
			end
			if isfield(structmd,'spcvx'), 
				md.stressbalance.spcvx=NaN*ones(md.mesh.numberofvertices,1);
				pos=find(~isnan(structmd.spcvx)); md.stressbalance.spcvx(pos)=structmd.spcvx(pos); 
			end
			if isfield(structmd,'spcvy'),
				md.stressbalance.spcvy=NaN*ones(md.mesh.numberofvertices,1);
				pos=find(~isnan(structmd.spcvy)); md.stressbalance.spcvy(pos)=structmd.spcvy(pos);     
			end
			if isfield(structmd,'spcvz'),
				md.stressbalance.spcvz=NaN*ones(md.mesh.numberofvertices,1);
				pos=find(~isnan(structmd.spcvz)); md.stressbalance.spcvz(pos)=structmd.spcvz(pos);     
			end
			if isfield(structmd,'pressureload'),
				if ~isempty(structmd.pressureload) & ismember(structmd.pressureload(end,end),[118 119 120]),
					pos=find(structmd.pressureload(:,end)==120); md.stressbalance.icefront(pos,end)=0;
					pos=find(structmd.pressureload(:,end)==118); md.stressbalance.icefront(pos,end)=1;
					pos=find(structmd.pressureload(:,end)==119); md.stressbalance.icefront(pos,end)=2;
				end
			end
			if isfield(structmd,'elements_type') & structmd.elements_type(end,end)>50,
				pos=find(structmd.elements_type==59); md.flowequation.element_equation(pos,end)=0;
				pos=find(structmd.elements_type==55); md.flowequation.element_equation(pos,end)=1;
				pos=find(structmd.elements_type==56); md.flowequation.element_equation(pos,end)=2;
				pos=find(structmd.elements_type==60); md.flowequation.element_equation(pos,end)=3;
				pos=find(structmd.elements_type==62); md.flowequation.element_equation(pos,end)=4;
				pos=find(structmd.elements_type==57); md.flowequation.element_equation(pos,end)=5;
				pos=find(structmd.elements_type==58); md.flowequation.element_equation(pos,end)=6;
				pos=find(structmd.elements_type==61); md.flowequation.element_equation(pos,end)=7;
			end
			if isfield(structmd,'vertices_type') & structmd.vertices_type(end,end)>50,
				pos=find(structmd.vertices_type==59); md.flowequation.vertex_equation(pos,end)=0;
				pos=find(structmd.vertices_type==55); md.flowequation.vertex_equation(pos,end)=1;
				pos=find(structmd.vertices_type==56); md.flowequation.vertex_equation(pos,end)=2;
				pos=find(structmd.vertices_type==60); md.flowequation.vertex_equation(pos,end)=3;
				pos=find(structmd.vertices_type==62); md.flowequation.vertex_equation(pos,end)=4;
				pos=find(structmd.vertices_type==57); md.flowequation.vertex_equation(pos,end)=5;
				pos=find(structmd.vertices_type==58); md.flowequation.vertex_equation(pos,end)=6;
				pos=find(structmd.vertices_type==61); md.flowequation.vertex_equation(pos,end)=7;
			end
			if isfield(structmd,'rheology_law') & isnumeric(structmd.rheology_law),
				if (structmd.rheology_law==272), md.materials.rheology_law='None';      end
				if (structmd.rheology_law==368), md.materials.rheology_law='Paterson';  end
				if (structmd.rheology_law==369), md.materials.rheology_law='Arrhenius'; end
			end
			if isfield(structmd,'groundingline_migration') & isnumeric(structmd.groundingline_migration),
				if (structmd.groundingline_migration==272), md.groundingline.migration='None';      end
				if (structmd.groundingline_migration==273), md.groundingline.migration='AggressiveMigration';  end
				if (structmd.groundingline_migration==274), md.groundingline.migration='SoftMigration'; end
			end
			if isfield(structmd,'control_type') & isnumeric(structmd.control_type),
				if (structmd.control_type==143), md.inversion.control_parameters={'FrictionCoefficient'}; end
				if (structmd.control_type==190), md.inversion.control_parameters={'RheologyBbar'}; end
				if (structmd.control_type==147), md.inversion.control_parameters={'Thickeningrate'}; end
			end
			if isfield(structmd,'cm_responses') & ismember(structmd.cm_responses(end,end),[165:170 383 388 389]),
				pos=find(structmd.cm_responses==166); md.inversion.cost_functions(pos)=101;
				pos=find(structmd.cm_responses==167); md.inversion.cost_functions(pos)=102;
				pos=find(structmd.cm_responses==168); md.inversion.cost_functions(pos)=103;
				pos=find(structmd.cm_responses==169); md.inversion.cost_functions(pos)=104;
				pos=find(structmd.cm_responses==170); md.inversion.cost_functions(pos)=105;
				pos=find(structmd.cm_responses==165); md.inversion.cost_functions(pos)=201;
				pos=find(structmd.cm_responses==389); md.inversion.cost_functions(pos)=501;
				pos=find(structmd.cm_responses==388); md.inversion.cost_functions(pos)=502;
				pos=find(structmd.cm_responses==382); md.inversion.cost_functions(pos)=503;
			end

			if isfield(structmd,'artificial_diffusivity') & structmd.artificial_diffusivity==2,
					md.thermal.stabilization=2;
					md.masstransport.stabilization=1;
					md.balancethickness.stabilization=1;
			end
			if isnumeric(md.masstransport.hydrostatic_adjustment)
				if md.masstransport.hydrostatic_adjustment==269,
					md.masstransport.hydrostatic_adjustment='Incremental';
				else
					md.masstransport.hydrostatic_adjustment='Absolute';
				end
			end

			%New fields
			if ~isfield(structmd,'upperelements') & isa(md.mesh,'mesh3dprisms')
				md.mesh.upperelements=transpose(1:md.mesh.numberofelements)+md.mesh.numberofelements2d;
				md.mesh.upperelements(end-md.mesh.numberofelements2d+1:end)=NaN;
			end
			if ~isfield(structmd,'lowerelements') & isa(md.mesh,'mesh3dprisms')
				md.mesh.lowerelements=transpose(1:md.mesh.numberofelements)-md.mesh.numberofelements2d;
				md.mesh.lowerelements(1:md.mesh.numberofelements2d)=NaN;
			end
			if ~isfield(structmd,'diagnostic_ref');
				md.stressbalance.referential=NaN*ones(md.mesh.numberofvertices,6);
			end
			if ~isfield(structmd,'loadingforce');
				md.stressbalance.loadingforce=0*ones(md.mesh.numberofvertices,3);
			end

			%2013 August 9
			if isfield(structmd,'prognostic') & isa(structmd.prognostic,'prognostic'),
				disp('Recovering old prognostic class');
				md.masstransport=masstransport(structmd.prognostic);
			end
			%2013 August 9
			if isfield(structmd,'diagnostic') & (isa(structmd.diagnostic,'diagnostic') || isa(structmd.diagnostic,'stressbalance')),
				disp('Recovering old diagnostic class');
				md.stressbalance=stressbalance(structmd.diagnostic);
			end
			%2014 January 9th
			if isfield(structmd,'surfaceforcings') & isa(md.smb,'surfaceforcings'),
				disp('Recovering old surfaceforcings class');
				mass_balance=structmd.surfaceforcings.mass_balance;
				md.smb=SMB();
				md.smb.mass_balance=mass_balance;
			end
			%2015 September 10
			if isfield(structmd,'surfaceforcings') & isa(structmd.surfaceforcings,'SMB'),
				disp('Recovering old SMB class');
				md.smb=SMBforcing(structmd.surfaceforcings);
			end
			if isfield(structmd,'surfaceforcings') & isa(structmd.surfaceforcings,'SMBhenning'),
				disp('Recovering old SMBhenning class');
				md.smb=SMBhenning(structmd.surfaceforcings);
			end
			if isfield(structmd,'slr') && ~isempty(structmd.slr)
				md.solidearth       = solidearth('Earth');
				disp('Recovering old slr class');
				if isfield(structmd.slr,'sealevel'),
					md.solidearth.sealevel=structmd.slr.sealevel;
				end
				md.solidearth.planetradius=structmd.slr.planetradius;
				md.solidearth.requested_outputs=structmd.slr.requested_outputs;
				md.solidearth.transitions=structmd.slr.transitions;

				md.solidearth.transitions=structmd.slr.transitions;
				md.solidearth.settings.reltol=structmd.slr.reltol;
				md.solidearth.settings.abstol=structmd.slr.abstol;
				md.solidearth.settings.maxiter=structmd.slr.maxiter;
				md.solidearth.settings.rigid=structmd.slr.rigid;
				md.solidearth.settings.elastic=structmd.slr.elastic;
				md.solidearth.settings.rotation=structmd.slr.rotation;
				md.solidearth.settings.runfrequency=structmd.slr.geodetic_run_frequency;
				md.solidearth.settings.computesealevelchange=structmd.slr.geodetic;
				md.solidearth.settings.degacc=structmd.slr.degacc;
				md.solidearth.settings.horiz=structmd.slr.horiz;
				md.solidearth.settings.ocean_area_scaling=structmd.slr.ocean_area_scaling;

				md.solidearth.surfaceload.icethicknesschange=structmd.slr.deltathickness;
				md.solidearth.surfaceload.waterheightchange=structmd.slr.hydro_rate;

				md.solidearth.lovenumbers.h=structmd.slr.love_h;
				md.solidearth.lovenumbers.k=structmd.slr.love_k;
				md.solidearth.lovenumbers.l=structmd.slr.love_l;
				md.solidearth.lovenumbers.th=structmd.slr.tide_love_h;
				md.solidearth.lovenumbers.tk=structmd.slr.tide_love_k;
				md.solidearth.lovenumbers.tk2secular=structmd.slr.fluid_love;

				md.solidearth.rotational.equatorialmoi=structmd.slr.equatorial_moi;
				md.solidearth.rotational.polarmoi=structmd.slr.polar_moi;
				md.solidearth.rotational.angularvelocity=structmd.slr.angular_velocity;
			end
		end% }}}
		function md = tetras(md,varargin) % {{{
			%TETRAS - split 3d prismatic mesh into 3 tetrahedrons
			%
			%   Usage:
			%      md=tetra(md)

			if ~isa(md.mesh,'mesh3dprisms')
				error('mesh is not a 3d prismatic mesh');
			end

			%Initialize tetra mesh
			md.mesh=mesh3dtetras(md.mesh);

			%Subdivision from Philipp Furnstahl (http://studierstube.icg.tugraz.at/thesis/fuernstahl_thesis.pdf)
			steiner  = 0;
			nbv      = md.mesh.numberofvertices;
			nbt      = 3*md.mesh.numberofelements;
			elements = zeros(nbt,4);
			for i=1:md.mesh.numberofelements
				v1=md.mesh.elements(i,1); v2=md.mesh.elements(i,2); v3=md.mesh.elements(i,3);
				v4=md.mesh.elements(i,4); v5=md.mesh.elements(i,5); v6=md.mesh.elements(i,6);
				if(min(v2,v4)<min(v1,v5) & min(v1,v6)<min(v3,v4) & min(v3,v5)<min(v2,v6)),
					steiner = steiner+1; nbv = nbv+1; nbt = nbt+5; v7 = nbv;
					md.mesh.x=[md.mesh.x; mean(md.mesh.x(md.mesh.elements(i,:)))];
					md.mesh.y=[md.mesh.y; mean(md.mesh.y(md.mesh.elements(i,:)))];
					md.mesh.z=[md.mesh.z; mean(md.mesh.z(md.mesh.elements(i,:)))];
					elements(3*(i-1)+1,:) = [v1 v2 v3 v7];
					elements(3*(i-1)+2,:) = [v1 v2 v4 v7];
					elements(3*(i-1)+3,:) = [v2 v4 v5 v7];
					elements(end+1,:) = [v2 v3 v5 v7];
					elements(end+1,:) = [v3 v5 v6 v7];
					elements(end+1,:) = [v1 v3 v6 v7];
					elements(end+1,:) = [v1 v4 v6 v7];
					elements(end+1,:) = [v4 v5 v6 v7];
				elseif(min(v2,v4)<min(v1,v5) & min(v1,v6)<min(v3,v4) & min(v3,v5)>min(v2,v6)),
					elements(3*(i-1)+1,:) = [v1 v2 v4 v6];
					elements(3*(i-1)+2,:) = [v2 v4 v5 v6];
					elements(3*(i-1)+3,:) = [v1 v2 v3 v6];
				elseif(min(v2,v4)<min(v1,v5) & min(v1,v6)>min(v3,v4) & min(v3,v5)<min(v2,v6)),
					elements(3*(i-1)+1,:) = [v1 v2 v3 v4];
					elements(3*(i-1)+2,:) = [v2 v3 v4 v5];
					elements(3*(i-1)+3,:) = [v3 v4 v5 v6];
				elseif(min(v2,v4)<min(v1,v5) & min(v1,v6)>min(v3,v4) & min(v3,v5)>min(v2,v6)),
					elements(3*(i-1)+1,:) = [v1 v2 v3 v4];
					elements(3*(i-1)+2,:) = [v2 v4 v5 v6];
					elements(3*(i-1)+3,:) = [v2 v3 v4 v6];
				elseif(min(v2,v4)>min(v1,v5) & min(v1,v6)<min(v3,v4) & min(v3,v5)<min(v2,v6)),
					elements(3*(i-1)+1,:) = [v1 v4 v5 v6];
					elements(3*(i-1)+2,:) = [v1 v2 v3 v5];
					elements(3*(i-1)+3,:) = [v1 v3 v5 v6];
				elseif(min(v2,v4)>min(v1,v5) & min(v1,v6)<min(v3,v4) & min(v3,v5)>min(v2,v6)),
					elements(3*(i-1)+1,:) = [v1 v4 v5 v6];
					elements(3*(i-1)+2,:) = [v1 v2 v5 v6];
					elements(3*(i-1)+3,:) = [v1 v2 v3 v6];
				elseif(min(v2,v4)>min(v1,v5) & min(v1,v6)>min(v3,v4) & min(v3,v5)<min(v2,v6)),
					elements(3*(i-1)+1,:) = [v1 v3 v4 v5];
					elements(3*(i-1)+2,:) = [v1 v2 v3 v5];
					elements(3*(i-1)+3,:) = [v3 v4 v5 v6];
				elseif(min(v2,v4)>min(v1,v5) & min(v1,v6)<min(v3,v4) & min(v3,v5)<min(v2,v6)),
					elements(3*(i-1)+1,:) = [v1 v5 v6 v4];
					elements(3*(i-1)+2,:) = [v1 v2 v3 v5];
					elements(3*(i-1)+3,:) = [v5 v6 v3 v1];
				elseif(min(v2,v4)>min(v1,v5) & min(v1,v6)>min(v3,v4) & min(v3,v5)>min(v2,v6)),
					steiner = steiner+1; nbv = nbv+1; nbt = nbt+5; v7 = nbv;
					md.mesh.x=[md.mesh.x; mean(md.mesh.x(md.mesh.elements(i,:)))];
					md.mesh.y=[md.mesh.y; mean(md.mesh.y(md.mesh.elements(i,:)))];
					md.mesh.z=[md.mesh.z; mean(md.mesh.z(md.mesh.elements(i,:)))];
					elements(3*(i-1)+1,:) = [v1 v2 v3 v7];
					elements(3*(i-1)+2,:) = [v1 v4 v5 v7];
					elements(3*(i-1)+3,:) = [v1 v2 v5 v7];
					elements(end+1,:) = [v2 v5 v6 v7];
					elements(end+1,:) = [v2 v3 v6 v7];
					elements(end+1,:) = [v3 v4 v6 v7];
					elements(end+1,:) = [v1 v3 v4 v7];
					elements(end+1,:) = [v4 v5 v6 v7];
				else
					error('Case not supported'); %not supposed to happen!
				end
				%Reorder elements to make sure they are direct
				for j=1:3
					element = elements(3*(i-1)+j,:);
					matrix = [md.mesh.x(element), md.mesh.y(element), md.mesh.z(element), ones(4,1)];
					if det(matrix)>0,
						elements(3*(i-1)+j,1)=element(2);
						elements(3*(i-1)+j,2)=element(1);
					end
				end
			end
			%%Split in 3 tetras
			%subelement1 = [1 2 3 5];
			%subelement2 = [4 6 5 1];
			%subelement3 = [5 6 3 1];
			%elements=[md.mesh.elements(:,subelement1);md.mesh.elements(:,subelement2);md.mesh.elements(:,subelement3)];
			if steiner==0,
				disp('No Steiner point required to split prismatic mesh into tets');
			else
				disp([num2str(steiner) ' Steiner points had to be included'])
				error('Steiner point not supported yet');
			end

			pos_elements = repmat([1:md.mesh.numberofelements]',3,1);

			md.mesh.elements=elements;
			md.mesh.numberofelements=size(elements,1);

			%p and q (same deal, except for element that are on the bedrock: )
			if ~isnan(md.friction.p),
				md.friction.p=md.friction.p(pos_elements);
				md.friction.q=md.friction.q(pos_elements);
			end

			%elementstype
			if ~isnan(md.flowequation.element_equation)
				oldelements_type=md.flowequation.element_equation;
				md.flowequation.element_equation=md.flowequation.element_equation(pos_elements);
			end

			%connectivity
			md.mesh.elementconnectivity=NaN;

			%materials
			if ~isnan(md.materials.rheology_n),
				md.materials.rheology_n=md.materials.rheology_n(pos_elements);
			end

			%increase connectivity if less than 25:
			if md.mesh.average_vertex_connectivity<=25,
				md.mesh.average_vertex_connectivity=100;
			end
		end % }}}
		function memory(self) % {{{

			disp(sprintf('\nMemory imprint:\n'));

			fields=properties('model');
			mem=0;

			for i=1:length(fields),
				field=self.(fields{i});
				s=whos('field'); 
				mem=mem+s.bytes/1e6;
				disp(sprintf('%19s: %6.2f Mb',fields{i},s.bytes/1e6));
			end
			disp(sprintf('%19s--%10s','--------------','--------------'));
			disp(sprintf('%19s: %g Mb','Total',mem));
		end
		% }}}
		function netcdf(self,filename) % {{{
			%NETCDF - save model as netcdf
			%
			%   Usage:
			%      netcdf(md,filename)
			%
			%   Example:
			%      netcdf(md,'model.nc');

			disp('Saving model as NetCDF');
			%1. Create NetCDF file
			ncid=netcdf.create(filename,'CLOBBER');
			netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Conventions','CF-1.4');
			netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Title',['ISSM model (' self.miscellaneous.name ')']);
			netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Author',getenv('USER'));
			netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Date',datestr(now));

			%Preallocate variable id, needed to write variables in netcdf file
			var_id=zeros(1000,1);%preallocate

			for step=1:2,
				counter=0;
				[var_id,counter]=structtonc(ncid,'md',self,0,var_id,counter,step);
				if step==1, netcdf.endDef(ncid); end
			end

			if counter>1000,
				warning(['preallocation of var_id need to be updated from ' num2str(1000) ' to ' num2str(counter)]);
			end

			netcdf.close(ncid)
		end % }}}
		function xylim(self) % {{{

			xlim([min(self.mesh.x) max(self.mesh.x)]);
			ylim([min(self.mesh.y) max(self.mesh.y)])
		end % }}}
		function md=upload(md) % {{{
			%the goal of this routine is to upload the model onto a server, and to empty it.
			%So first, save the model with a unique name and upload the file to the server: 
			random_part=fix(rand(1)*10000);
			id=[md.miscellaneous.name '-' regexprep(datestr(now),'[^\w'']','') '-' num2str(random_part)  '-' getenv('USER') '-' oshostname() '.upload']; 
			save('id','md');

			%Now, upload the file: 
			issmscpout(md.settings.upload_server,md.settings.upload_path,md.settings.upload_login,md.settings.upload_port,{id},1);

			%Now, empty this model of everything except settings, and record name of file we just uploaded!
			settings_back=md.settings;
			md=model();
			md.settings=settings_back;
			md.settings.upload_filename=id;

			%get locally rid of file that was uploaded
			delete(id);

		end % }}}
		function md=download(md) % {{{

			%the goal of this routine is to download the internals of the current model from a server, because 
			%this model is empty, except for the settings which tell us where to go and find this model!

			%Download the file: 
			issmscpin(md.settings.upload_server, md.settings.upload_login, md.settings.upload_port, md.settings.upload_path, {md.settings.upload_filename});

			name=md.settings.upload_filename;

			%Now, load this model: 
			md=loadmodel(md.settings.upload_filename);

			%get locally rid of file that was downloaded
			delete(name);

		end % }}}
		function saveasstruct(md,filename) % {{{

			%Get all model fields
			mdfields=properties('model');

			disp('Converting all model fields to struct...');
			warning off MATLAB:structOnObject
			for i=1:length(mdfields),

				%convert md field to struct
				field=mdfields{i};
				field_struct = struct(md.(field));

				%Check that there is no class remaining in the field
				subfields = fields(field_struct);
				for i=1:numel(subfields)
					if isobject(field_struct.(subfields{i}))
						disp(['skipping ' subfields{i} ' because it is an object'])
						field_struct = rmfield(field_struct, subfields{i});
					end
				end
				md.(field) = field_struct;
			end

			disp('Converting model to struct...');
			md=struct(md);

			disp(['Saving as ' filename '...']);
			warning on MATLAB:structOnObject
			save(filename,'md','-v7.3')

		end % }}}
		function savemodeljs(md,modelname,websiteroot,varargin) % {{{

			%the goal of this routine is to save the model as a javascript array that can be included in any html 
			%file: 

			options=pairoptions(varargin{:});
			optimization=getfieldvalue(options,'optimize',0);


			%disp: 
			disp(['saving model ''' modelname ''' in file ' websiteroot '/js/' modelname '.js']);

			%open file for writing and declare the model:
			fid=fopen([websiteroot '/js/' modelname '.js'],'w');
			fprintf(fid,'var %s=new model();\n',modelname);

			%now go through all the classes and fwrite all the corresponding fields: 

			fields=properties('model');
			for i=1:length(fields),
				field=fields{i};

				%Some properties do not need to be saved
				if ismember(field,{'results','cluster'}),
					continue;
				end

				%some optimization: 
				if optimization==1,
					%optimize for plotting only:
					if ~ismember(field,{'geometry','mesh','mask'}),
						continue;
					end
				end

				%Check that current field is an object
				if ~isobject(md.(field))
					error(['field ''' char(field) ''' is not an object']);
				end

				if ~ismethod(md.(field),'savemodeljs')
					disp(['Note: Method ''savemodeljs'' not yet implemented for class ' field])
				else
					%savemodeljs for current object
					%disp(['javascript saving ' field '...']);
					savemodeljs(md.(field),fid,modelname);
				end
			end

			%done, close file:
			fclose(fid);
		end
	end
end
