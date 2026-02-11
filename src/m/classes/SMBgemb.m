%SMBgemb Class definition.
%   This is the class that hosts all the inputs for the Alberta Glacier Surface Mass Balance Model
%   Alex Gardner, University of Alberta.
%
%   Usage:
%      SMBgemb=SMBgemb(md.mesh);

classdef SMBgemb
	properties (SetAccess=public)
		% {{{
		%each one of these properties is a transient forcing to the GEMB model, loaded from meteorological data derived
		%from an automatic weather stations (AWS). Each property is therefore a matrix, of size (numberofvertices x number
		%of time steps. )

		%solution choices
		isgraingrowth       = 0;
		isalbedo            = 0;
		isshortwave         = 0;
		isthermal           = 0;
		isaccumulation      = 0;
		ismelt              = 0;
		isdensification     = 0;
		isturbulentflux     = 0;
		isconstrainsurfaceT = 0;
		isdeltaLWup         = 0;
		ismappedforcing     = 0;
		isprecipforcingremapped = 0;
		iscompressedforcing = 0;

		%inputs:
		Ta    = NaN; %2 m air temperature, in Kelvin
		V     = NaN; %wind speed (m/s-1)
		dswrf = NaN; %downward shortwave radiation flux [W/m^2]
		dlwrf = NaN; %downward longwave radiation flux [W/m^2]
		P     = NaN; %precipitation [mm w.e. / m^2]
		eAir  = NaN; %screen level vapor pressure [Pa]
		pAir  = NaN; %surface pressure [Pa]

		Tmean = NaN; %mean annual temperature [K]
		Vmean = NaN; %mean annual wind velocity [m s-1]
		C     = NaN; %mean annual snow accumulation [kg m-2 yr-1]
		Tz    = NaN; %height above ground at which temperature (T) was sampled [m]
		Vz    = NaN; %height above ground at which wind (V) was sampled [m]

		%optional inputs:
		aValue = NaN; %Albedo forcing at every element.  Used only if aIdx == 0, or density exceeds adThresh
		teValue = NaN; %Outward longwave radiation thermal emissivity forcing at every element (default in code is 1). 
		               %Used only if eIdx== 0, or effective grain radius exceeds teThresh
		dulwrfValue = NaN; %Delta with which to perturb the long wave radiation upwards. Use if isdeltaLWup is true.  
		mappedforcingpoint = NaN; %Mapping of which forcing point will map to each mesh element (integer). Of size number of elements.
		                        %Use if ismappedforcing is true.
		mappedforcingelevation = NaN; %The elevation of each mapped forcing location (m above sea level). Of size number
		                        %of forcing points. Use if ismappedforcing is true.
		lapseTaValue = NaN; %Temperature lapse rate if forcing has different grid and should be remapped. Use if ismappedforcing is true.
								  % (Default value is -0.006 K m-1., vector of mapping points)
		lapsedlwrfValue = NaN; %Longwave down lapse rate if forcing has different grid and should be remapped. Use if ismappedforcing is true.
		                    % Where set to 0, dlwrf will scale with a constant effective atmospheric emissivity. 
								  % (Default value is -0.032 W m-2 m-1., vector of mapping points)

		% Initialization of snow properties
		Dzini = NaN; %cell depth (m)
		Dini = NaN; %snow density (kg m-3)
		Reini = NaN; %effective grain size (mm)
		Gdnini = NaN; %grain dendricity (0-1)
		Gspini = NaN; %grain sphericity (0-1)
		ECini = NaN; %evaporation/condensation (kg m-2)
		Wini = NaN; %Water content (kg m-2)
		Aini = NaN; %albedo (0-1)
		Adiffini = NaN; %albedo, diffusive radiation (0-1)
		Tini = NaN; %snow temperature (K)
		Sizeini = NaN; %Number of layers

		%settings:
		aIdx   = NaN; %method for calculating albedo and subsurface absorption (default is 1)
		% 0: direct input from aValue parameter, no use of adThresh
		% 1: effective grain radius [Gardner & Sharp, 2009]
		% 2: effective grain radius [Brun et al., 1992; LeFebre et al., 2003]], with swIdx=1, SW penetration follows grain size in 3 spectral bands (Brun et al., 1992)
		% 3: density and cloud amount [Greuell & Konzelmann, 1994]
		% 4: exponential time decay & wetness [Bougamont & Bamber, 2005]

		eIdx   = NaN; %method for calculating emissivity (default is 1)
		% 0: direct input from teValue parameter, no use of teThresh
		% 1: default value of 1, in areas with grain radius below teThresh
		% 2: default value of 1, in areas with grain radius below teThresh and areas of dry snow (not bare ice or wet) at the surface

		tcIdx   = NaN; %method for calculating thermal conductivity (default is 1)
		% 1: after Sturm et al, 1997
		% 2: after Calonne et al., 2011

		swIdx  = NaN; %apply all SW to top grid cell (0) or allow SW to penetrate surface (1) (default 0, if swIdx=1 and aIdx=2, function of effective radius (Brun et al., 1992) or else dependent on snow density (taken from Bassford, 2002))

		denIdx = NaN; %densification model to use (default is 2):
		% 1 = emperical model of Herron and Langway (1980)
		% 2 = semi-emperical model of Anthern et al. (2010)
		% 3 = DO NOT USE: physical model from Appendix B of Anthern et al. (2010)
		% 4 = DO NOT USE: emperical model of Li and Zwally (2004)
		% 5 = DO NOT USE: modified emperical model (4) by Helsen et al. (2008)
		% 6 = Antarctica semi-emperical model of Ligtenberg et al. (2011)
		% 7 = Greenland semi-emperical model of Kuipers Munneke et al. (2015)

		dsnowIdx = NaN; %model for fresh snow accumulation density (default is 1):
		% 0 = Original GEMB value, 150 kg/m^3
		% 1 = Antarctica value of fresh snow density, 350 kg/m^3
		% 2 = Greenland value of fresh snow density, 315 kg/m^3, Fausto et al. (2018)
		% 3 = Antarctica model of Kaspers et al. (2004)
		% 4 = Greenland model of Kuipers Munneke et al. (2015)

		zTop  = NaN; % depth over which grid length is constant at the top of the snopack (default 10) [m]
		dzTop = NaN; % initial top vertical grid spacing (default .05) [m]
		dzMin = NaN; % initial min vertical allowable grid spacing (default dzMin/2) [m]

		zY    = NaN; % stretch grid cells bellow top_z by a [top_dz * y ^ (cells bellow top_z)]
		zMax = NaN; %initial max model depth (default is min(thickness,250)) [m]
		zMin = NaN; %initial min model depth (default is min(thickness,130)) [m]
		outputFreq = NaN; %output frequency in days (default is monthly, 30)

		%specific albedo parameters:
		%Method 1
		dswdiffrf = NaN; %downward diffusive shortwave radiation flux [W/m^2]
		szaValue = NaN; %Solar Zenith Angle [degree]
		cotValue = NaN; %Cloud Optical Thickness
		ccsnowValue = NaN; %concentration of light absorbing carbon for snow [ppm1]
		cciceValue = NaN; %concentration of light absorbing carbon for ice [ppm1]
		%Method 1 and 2:
		aSnow = NaN; % new snow albedo (0.64 - 0.89)
		aIce  = NaN; % range 0.27-0.58 for old snow
		%Method 3: Radiation Correction Factors -> only used for met station data and Greuell & Konzelmann, 1994 albedo
		cldFrac = NaN; % average cloud amount
		%Method 4: additonal tuning parameters albedo as a funtion of age and water content (Bougamont et al., 2005)
		t0wet = NaN; % time scale for wet snow (15-21.9)
		t0dry = NaN; % warm snow timescale (30)
		K     = NaN; % time scale temperature coef. (7)
		adThresh = NaN; %Apply aIdx method to all areas with densities below this value,
		%or else apply direct input value from aValue, allowing albedo to be altered.
		%Default value is rho water (1023 kg m-3).
		teThresh = NaN; %Apply eIdx method to all areas with grain radii above this value (mm),
		%or else apply direct input value from teValue, allowing emissivity to be altered.
		%Default value is a effective grain radius of 10 mm.

		%densities:
		InitDensityScaling= NaN; %initial scaling factor multiplying the density of ice, which describes the density of the snowpack.

		%thermal:
		ThermoDeltaTScaling= NaN; %scaling factor to multiply the thermal diffusion timestep (delta t)

		steps_per_step = 1;
		averaging = 0;
		requested_outputs      = {};

		%Several fields are missing from the standard GEMB model, which are 
		%captured intrinsically by ISSM.
		%dateN: that's the last row of the above fields.
		%dt:    included in dateN. Not an input.
		%elev:  this is taken from the ISSM surface itself.

	end % }}}
	methods
		function self = SMBgemb(varargin) % {{{
			switch nargin
				case 1
					mesh=varargin{1};
					self=setdefaultparameters(self,mesh);
				otherwise
					error('constructor not supported: need mesh to set defaults');
			end
		end % }}}
		function disp(self) % {{{

			disp(sprintf('   surface forcings for SMB GEMB model :'));

			fielddisplay(self,'isgraingrowth','run grain growth module (default true)');
			fielddisplay(self,'isalbedo','run albedo module (default true)');
			fielddisplay(self,'isshortwave','run short wave module (default true)');
			fielddisplay(self,'isthermal','run thermal module (default true)');
			fielddisplay(self,'isaccumulation','run accumulation module (default true)');
			fielddisplay(self,'ismelt','run melting  module (default true)');
			fielddisplay(self,'isdensification','run densification module (default true)');
			fielddisplay(self,'isturbulentflux','run turbulant heat fluxes module (default true)');
			fielddisplay(self,'isconstrainsurfaceT','constrain surface temperatures to air temperature, turn off EC and surface flux contribution to surface temperature change (default false)');
			fielddisplay(self,'isdeltaLWup','set to true to invoke a bias in the long wave upward spatially, specified by dulwrfValue (default false)'); 
			fielddisplay(self,'ismappedforcing','set to true if forcing grid does not match model mesh, mapping specified by mappedforcingpoint (default false)');
			fielddisplay(self,'isprecipforcingremapped','set to true if ismappedforcing is true and precip should be downscaled from native grid (Default value is true)');
			fielddisplay(self,'iscompressedforcing','set to true to compress the input matrices when writing to binary (default false)');
			fielddisplay(self,'Ta','2 m air temperature, in Kelvin');
			fielddisplay(self,'V','wind speed (m s-1)');
			fielddisplay(self,'dswrf','downward shortwave radiation flux [W/m^2]');
			fielddisplay(self,'dswdiffrf','downward diffusive portion of shortwave radiation flux (default to 0) [W/m^2]');
			fielddisplay(self,'dlwrf','downward longwave radiation flux [W/m^2]');
			fielddisplay(self,'P','precipitation [mm w.e. / m^2]');
			fielddisplay(self,'eAir','screen level vapor pressure [Pa]');
			fielddisplay(self,'pAir','surface pressure [Pa]');
			fielddisplay(self,'Tmean','mean annual temperature [K]');
			fielddisplay(self,'C','mean annual snow accumulation [kg m-2 yr-1]');
			fielddisplay(self,'Vmean','mean annual wind speed [m s-1] (default 10 m/s)');
			fielddisplay(self,'Tz','height above ground at which temperature (T) was sampled [m]');
			fielddisplay(self,'Vz','height above ground at which wind (V) was sampled [m]');
			fielddisplay(self,'zTop','depth over which grid length is constant at the top of the snopack (default 10) [m]');
			fielddisplay(self,'dzTop','initial top vertical grid spacing (default .05) [m] ');
			fielddisplay(self,'dzMin','initial min vertical allowable grid spacing (default dzMin/2) [m] ');
			fielddisplay(self,'zMax','initial max model depth (default is min(thickness,250)) [m]');
			fielddisplay(self,'zMin','initial min model depth (default is min(thickness,130)) [m]');
			fielddisplay(self,'zY','stretch grid cells bellow top_z by a [top_dz * y ^ (cells bellow top_z)]');
			fielddisplay(self,'InitDensityScaling',{'initial scaling factor multiplying the density of ice','which describes the density of the snowpack.'});
			fielddisplay(self,'ThermoDeltaTScaling',{'scaling factor to multiply the thermal diffusion timestep (delta t)'});
			fielddisplay(self,'outputFreq','output frequency in days (default is monthly, 30)');
			fielddisplay(self,'adThresh','Apply aIdx method to all areas with densities below this value, or else apply direct input value from aValue, allowing albedo to be altered.');
			fielddisplay(self,'aIdx',{'method for calculating albedo and subsurface absorption (default is 1)',...
				'0: direct input from aValue parameter',...
				'1: effective grain radius [Gardner & Sharp, 2009]',...
				'2: effective grain radius [Brun et al., 1992; LeFebre et al., 2003], with swIdx=1, SW penetration follows grain size in 3 spectral bands (Brun et al., 1992)',...
				'3: density and cloud amount [Greuell & Konzelmann, 1994]',...
				'4: exponential time decay & wetness [Bougamont & Bamber, 2005]'})

			fielddisplay(self,'dulwrfValue','Specified bias to be applied to the outward long wave radiation at every element (W/m-2, +upward)');
			fielddisplay(self,'teValue','Outward longwave radiation thermal emissivity forcing at every element (default in code is 1)');
			fielddisplay(self,'teThresh',{'Apply eIdx method to all areas with effective grain radius above this value (mm),','or else apply direct input value from teValue, allowing emissivity to be altered.'});
			fielddisplay(self,'eIdx',{'method for calculating emissivity (default is 1)',...
				'0: direct input from teValue parameter, no use of teThresh',...
				'1: default value of 1, in areas with grain radius below teThresh',...
				'2: default value of 1, in areas with grain radius below teThresh and areas of dry snow (not bare ice or wet) at the surface'});

			fielddisplay(self,'tcIdx',{'method for calculating thermal conductivity (default is 1)',...
				'1: after Sturm et al, 1997',...
				'2: after Calonne et al., 2011'});

			fielddisplay(self,'mappedforcingpoint','Mapping of which forcing point will map to each mesh element for ismappedforcing option (integer). Size number of elements.');
			fielddisplay(self,'mappedforcingelevation','The elevation of each mapped forcing location (m above sea level) for ismappedforcing option. Size number of forcing points.');
			fielddisplay(self,'lapseTaValue','Temperature lapse rate of each mapped forcing location, if forcing has different grid and should be remapped for ismappedforcing option. (Default value is -0.006 K m-1, vector of mapping points)');
			fielddisplay(self,'lapsedlwrfValue','Longwave down lapse rate of each mapped forcing location, if forcing has different grid and should be remapped for ismappedforcing option. Where set to 0, dlwrf will scale with a constant effective atmospheric emissivity. (Default value is -0.032 W m-2 m-1, vector of mapping points)');

			%snow properties init
			fielddisplay(self,'Dzini','Initial cell depth when restart [m]');
			fielddisplay(self,'Dini','Initial snow density when restart [kg m-3]');
			fielddisplay(self,'Reini','Initial grain size when restart [mm]');
			fielddisplay(self,'Gdnini','Initial grain dendricity when restart [-]');
			fielddisplay(self,'Gspini','Initial grain sphericity when restart [-]');
			fielddisplay(self,'ECini','Initial evaporation/condensation when restart [kg m-2]');
			fielddisplay(self,'Wini','Initial snow water content when restart [kg m-2]');
			fielddisplay(self,'Aini','Initial albedo when restart [-]');
			fielddisplay(self,'Adiffini','Initial diffusive radiation albedo when restart (default to 1) [-]');
			fielddisplay(self,'Tini','Initial snow temperature when restart [K]');
			fielddisplay(self,'Sizeini','Initial number of layers when restart [-]');

			%additional albedo parameters:
			fielddisplay(self,'aValue','Albedo forcing at every element');
			switch self.aIdx
				case {1 2}
					fielddisplay(self,'aSnow','new snow albedo (0.64 - 0.89)');
					fielddisplay(self,'aIce','albedo of ice (0.27-0.58)');
					if self.aIdx==1
						fielddisplay(self,'szaValue','Solar Zenith Angle [degree]');
						fielddisplay(self,'cotValue','Cloud Optical Thickness');
						fielddisplay(self,'ccsnowValue','concentration of light absorbing carbon for snow [ppm1]');
						fielddisplay(self,'cciceValue','concentration of light absorbing carbon for ice [ppm1]');
					end
				case 3
					fielddisplay(self,'cldFrac','average cloud amount');
				case 4
					fielddisplay(self,'t0wet','time scale for wet snow (15-21.9) [d]');
					fielddisplay(self,'t0dry','warm snow timescale (30) [d]');
					fielddisplay(self,'K','time scale temperature coef. (7) [d]');
			end

			fielddisplay(self,'swIdx','apply all SW to top grid cell (0) or allow SW to penetrate surface (1) [default 0, if swIdx=1 and aIdx=2 function of effective radius (Brun et al., 1992) or else dependent on snow density (taken from Bassford, 2002)]');
			fielddisplay(self,'denIdx',{'densification model to use (default is 2):',...
				'1 = emperical model of Herron and Langway (1980)',...
				'2 = semi-emperical model of Anthern et al. (2010)',...
				'3 = DO NOT USE: physical model from Appendix B of Anthern et al. (2010)',...
				'4 = DO NOT USE: emperical model of Li and Zwally (2004)',...
				'5 = DO NOT USE: modified emperical model (4) by Helsen et al. (2008)',...
				'6 = Antarctica semi-emperical model of Ligtenberg et al. (2011)',...
				'7 = Greenland semi-emperical model of Kuipers Munneke et al. (2015)'});
			fielddisplay(self,'dsnowIdx',{'model for fresh snow accumulation density (default is 1):',...
				'0 = Original GEMB value, 150 kg/m^3',...
				'1 = Antarctica value of fresh snow density, 350 kg/m^3',...
				'2 = Greenland value of fresh snow density, 315 kg/m^3, Fausto et al. (2018)',...
				'3 = Antarctica model of Kaspers et al. (2004), Make sure to set Vmean accurately',...
				'4 = Greenland model of Kuipers Munneke et al. (2015)'});

			fielddisplay(self, 'steps_per_step', 'number of smb steps per time step');
			fielddisplay(self, 'averaging', 'averaging methods from short to long steps');
			disp(sprintf('%51s  0: Arithmetic (default)',' '));
			disp(sprintf('%51s  1: Geometric',' '));
			disp(sprintf('%51s  2: Harmonic',' '));
			fielddisplay(self,'requested_outputs','additional outputs requested');


		end % }}}
		function self = extrude(self,md) % {{{

			if size(self.Ta,1)==md.mesh.numberofelements | size(self.Ta,1)==md.mesh.numberofelements+1
				self.Ta=project3d(md,'vector',self.Ta,'type','element');
				self.V=project3d(md,'vector',self.V,'type','element');
				self.dswrf=project3d(md,'vector',self.dswrf,'type','element');
				self.dlwrf=project3d(md,'vector',self.dlwrf,'type','element');
				self.P=project3d(md,'vector',self.P,'type','element');
				self.eAir=project3d(md,'vector',self.eAir,'type','element');
				self.pAir=project3d(md,'vector',self.pAir,'type','element');
				self.C=project3d(md,'vector',self.C,'type','element');
				self.Tmean=project3d(md,'vector',self.Tmean,'type','element');
			end

			if ~isnan(self.Dzini)
				self.Dzini=project3d(md,'vector',self.Dzini,'type','element');
			end
			if ~isnan(self.Dini)
				self.Dini=project3d(md,'vector',self.Dini,'type','element');
			end
			if ~isnan(self.Reini)
				self.Reini=project3d(md,'vector',self.Reini,'type','element');
			end
			if ~isnan(self.Gdnini)
				self.Gdnini=project3d(md,'vector',self.Gdnini,'type','element');
			end
			if ~isnan(self.Gspini)
				self.Gspini=project3d(md,'vector',self.Gspini,'type','element');
			end
			if ~isnan(self.ECini)
				self.ECini=project3d(md,'vector',self.ECini,'type','element');
			end
			if ~isnan(self.Wini)
				self.Wini=project3d(md,'vector',self.Wini,'type','element');
			end
			if ~isnan(self.Aini)
				self.Aini=project3d(md,'vector',self.Aini,'type','element');
			end
			if ~isnan(self.Adiffini)
				self.Adiffini=project3d(md,'vector',self.Adiffini,'type','element');
			end
			if ~isnan(self.Tini)
				self.Tini=project3d(md,'vector',self.Tini,'type','element');
			end

			if ~isnan(self.dswdiffrf)
				self.dswdiffrf=project3d(md,'vector',self.dswdiffrf,'type','element');
			end
			if ~isnan(self.szaValue)
				self.szaValue=project3d(md,'vector',self.szaValue,'type','element');
			end
			if ~isnan(self.cotValue)
				self.cotValue=project3d(md,'vector',self.cotValue,'type','element');
			end
			if ~isnan(self.ccsnowValue)
				self.ccsnowValue=project3d(md,'vector',self.ccsnowValue,'type','element');
			end
			if ~isnan(self.cciceValue)
				self.cciceValue=project3d(md,'vector',self.cciceValue,'type','element');
			end
			if ~isnan(self.aValue)
				self.aValue=project3d(md,'vector',self.aValue,'type','element');
			end
			if ~isnan(self.teValue)
				self.teValue=project3d(md,'vector',self.teValue,'type','element');
			end
			if ~isnan(self.mappedforcingpoint)
				self.mappedforcingpoint=project3d(md,'vector',self.mappedforcingpoint,'type','element');
			end

		end % }}}
		function list = defaultoutputs(self,md) % {{{
			list = {'SmbMassBalance','SmbAccumulatedMassBalance'};
		end % }}}
		function self = setdefaultparameters(self,mesh) % {{{

			self.isgraingrowth=1;
			self.isalbedo=1;
			self.isshortwave=1;
			self.isthermal=1;
			self.isaccumulation=1;
			self.ismelt=1;
			self.isdensification=1;
			self.isturbulentflux=1;
			self.isconstrainsurfaceT=0;
			self.isdeltaLWup=0;
			self.ismappedforcing=0;
			self.isprecipforcingremapped=1;
			self.iscompressedforcing=0;

			self.aIdx = 1;
			self.eIdx = 1;
			self.tcIdx = 1;
			self.swIdx = 0;
			self.denIdx = 2;
			self.dsnowIdx = 1;
			self.zTop=10*ones(mesh.numberofelements,1);
			self.dzTop = .05* ones (mesh.numberofelements,1);
			self.dzMin = self.dzTop/2;
			self.InitDensityScaling = 1.0;
			self.ThermoDeltaTScaling = 1/11;

			self.Vmean=10.0*ones(mesh.numberofelements,1);

			self.zMax=250*ones(mesh.numberofelements,1);
			self.zMin=130*ones(mesh.numberofelements,1);
			self.zY = 1.025*ones(mesh.numberofelements,1);
			self.outputFreq = 30;

			%additional albedo parameters
			self.aSnow = 0.85;
			self.aIce = 0.48;
			self.cldFrac = 0.1;
			self.t0wet = 15;
			self.t0dry = 30;
			self.K = 7;
			self.adThresh = 1023;
			self.teThresh = 10;

			self.teValue = ones(mesh.numberofelements,1);
			self.aValue = self.aSnow*ones(mesh.numberofelements,1);
			self.dulwrfValue = zeros(mesh.numberofelements,1);

			self.lapseTaValue = -0.006;
			self.lapsedlwrfValue = -0.032; 

			self.dswdiffrf=0.0*ones(mesh.numberofelements,1);
			self.szaValue=0.0*ones(mesh.numberofelements,1);
			self.cotValue=0.0*ones(mesh.numberofelements,1);
			self.ccsnowValue=0.0*ones(mesh.numberofelements,1);
			self.cciceValue=0.0*ones(mesh.numberofelements,1);

			self.Dzini=0.05*ones(mesh.numberofelements,2);
			self.Dini=910.0*ones(mesh.numberofelements,2);
			self.Reini=2.5*ones(mesh.numberofelements,2);
			self.Gdnini=0.0*ones(mesh.numberofelements,2);
			self.Gspini=0.0*ones(mesh.numberofelements,2);
			self.ECini=0.0*ones(mesh.numberofelements,1);
			self.Wini=0.0*ones(mesh.numberofelements,2);
			self.Aini=self.aSnow*ones(mesh.numberofelements,2);
			self.Adiffini=ones(mesh.numberofelements,2);
			self.Tini=273.15*ones(mesh.numberofelements,2);
			% 		/!\ Default value of Tini must be equal to Tmean but don't know Tmean yet (computed when atmospheric forcings are interpolated on mesh).
			% 		If initialization without restart, this value will be overwritten when snow parameters are retrieved in Element.cpp
			self.Sizeini=2*ones(mesh.numberofelements,1);

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			md = checkfield(md,'fieldname','smb.isgraingrowth','values',[0 1]);
			md = checkfield(md,'fieldname','smb.isalbedo','values',[0 1]);
			md = checkfield(md,'fieldname','smb.isshortwave','values',[0 1]);
			md = checkfield(md,'fieldname','smb.isthermal','values',[0 1]);
			md = checkfield(md,'fieldname','smb.isaccumulation','values',[0 1]);
			md = checkfield(md,'fieldname','smb.ismelt','values',[0 1]);
			md = checkfield(md,'fieldname','smb.isdensification','values',[0 1]);
			md = checkfield(md,'fieldname','smb.isturbulentflux','values',[0 1]);
			md = checkfield(md,'fieldname','smb.isconstrainsurfaceT','values',[0 1]);
			md = checkfield(md,'fieldname','smb.isdeltaLWup','values',[0 1]);
			md = checkfield(md,'fieldname','smb.ismappedforcing','values',[0 1]);
			md = checkfield(md,'fieldname','smb.isprecipforcingremapped','values',[0 1]);
			md = checkfield(md,'fieldname','smb.iscompressedforcing','values',[0 1]);

			sizeta=size(self.Ta);
			md = checkfield(md,'fieldname','smb.Ta','mappedtimeseries',1,'NaN',1,'Inf',1,'>',273-100,'<',273+100); %-100/100 celsius min/max value
			md = checkfield(md,'fieldname','smb.V','mappedtimeseries',1,'NaN',1,'Inf',1,'>=',0,'<',45,'size',sizeta); %max 500 km/h
			md = checkfield(md,'fieldname','smb.dswrf','mappedtimeseries',1,'NaN',1,'Inf',1,'>=',0,'<=',1400,'size',sizeta);
			md = checkfield(md,'fieldname','smb.dswdiffrf','mappedtimeseries',1,'NaN',1,'Inf',1,'>=',0,'<=',1400);
			md = checkfield(md,'fieldname','smb.dlwrf','mappedtimeseries',1,'NaN',1,'Inf',1,'>=',0,'size',sizeta);
			md = checkfield(md,'fieldname','smb.P','mappedtimeseries',1,'NaN',1,'Inf',1,'>=',0,'<=',200,'size',sizeta);
			md = checkfield(md,'fieldname','smb.eAir','mappedtimeseries',1,'NaN',1,'Inf',1,'size',sizeta);

			md = checkfield(md,'fieldname','smb.Tmean','size',[sizeta(1)-1 1],'NaN',1,'Inf',1,'>',273-100,'<',273+100); %-100/100 celsius min/max value
			md = checkfield(md,'fieldname','smb.C','size',[sizeta(1)-1 1],'NaN',1,'Inf',1,'>',0);
			md = checkfield(md,'fieldname','smb.Vmean','size',[sizeta(1)-1 1],'NaN',1,'Inf',1,'>=',0);
			md = checkfield(md,'fieldname','smb.Tz','size',[sizeta(1)-1 1],'NaN',1,'Inf',1,'>=',0,'<=',5000);
			md = checkfield(md,'fieldname','smb.Vz','size',[sizeta(1)-1 1],'NaN',1,'Inf',1,'>=',0,'<=',5000);

			md = checkfield(md,'fieldname','smb.teValue','timeseries',1,'NaN',1,'Inf',1,'>=',0,'<=',1);
			md = checkfield(md,'fieldname','smb.dulwrfValue','timeseries',1,'NaN',1,'Inf',1);

			if (self.ismappedforcing)
				md = checkfield(md,'fieldname','smb.mappedforcingpoint','size',[md.mesh.numberofelements 1],'NaN',1,'Inf',1,'>',0,'<=',sizeta(1)-1);
				md = checkfield(md,'fieldname','smb.mappedforcingelevation','size',[sizeta(1)-1 1],'NaN',1,'Inf',1);
				if prod(size(self.lapseTaValue))==1
					disp('WARNING:smb.lapseTaValue is now a vector of mapped elements. Set to md.smb.lapseTaValue*ones(size(md.smb.mappedforcingelevation)).');
				end
				if prod(size(self.lapsedlwrfValue))==1
					disp('WARNING:smb.lapsedlwrfValue is now a vector of mapped elements. Set to md.smb.lapsedlwrfValue*ones(size(md.smb.mappedforcingelevation)).');
				end
				md = checkfield(md,'fieldname','smb.lapseTaValue','size',[sizeta(1)-1 1],'NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','smb.lapsedlwrfValue','size',[sizeta(1)-1 1], 'NaN',1,'Inf',1);
			end

			md = checkfield(md,'fieldname','smb.aIdx','NaN',1,'Inf',1,'values',[0,1,2,3,4]);
			md = checkfield(md,'fieldname','smb.eIdx','NaN',1,'Inf',1,'values',[0,1,2]);
			md = checkfield(md,'fieldname','smb.tcIdx','NaN',1,'Inf',1,'values',[1,2]);
			md = checkfield(md,'fieldname','smb.swIdx','NaN',1,'Inf',1,'values',[0,1]);
			md = checkfield(md,'fieldname','smb.denIdx','NaN',1,'Inf',1,'values',[1,2,3,4,5,6,7]);
			md = checkfield(md,'fieldname','smb.dsnowIdx','NaN',1,'Inf',1,'values',[0,1,2,3,4]);

			md = checkfield(md,'fieldname','smb.zTop','NaN',1,'Inf',1,'>=',0);
			md = checkfield(md,'fieldname','smb.dzTop','NaN',1,'Inf',1,'>',0);
			md = checkfield(md,'fieldname','smb.dzMin','NaN',1,'Inf',1,'>',0);
			md = checkfield(md,'fieldname','smb.zY','NaN',1,'Inf',1,'>=',1);
			md = checkfield(md,'fieldname','smb.outputFreq','NaN',1,'Inf',1,'>',0,'<',10*365); %10 years max
			md = checkfield(md,'fieldname','smb.InitDensityScaling','NaN',1,'Inf',1,'>=',0,'<=',1);
			md = checkfield(md,'fieldname','smb.ThermoDeltaTScaling','NaN',1,'Inf',1,'>=',0,'<=',1);
			md = checkfield(md,'fieldname','smb.adThresh','NaN',1,'Inf',1,'>=',0);
			md = checkfield(md,'fieldname','smb.teThresh','NaN',1,'Inf',1,'>=',0);

			md = checkfield(md,'fieldname','smb.aValue','timeseries',1,'NaN',1,'Inf',1,'>=',0,'<=',1);
			switch self.aIdx,
				case {1 2}
					md = checkfield(md,'fieldname','smb.aSnow','NaN',1,'Inf',1,'>=',.64,'<=',.89);
					md = checkfield(md,'fieldname','smb.aIce','NaN',1,'Inf',1,'>=',.27,'<=',.58);
					if self.aIdx==1
						md = checkfield(md,'fieldname','smb.szaValue','timeseries',1,'NaN',1,'Inf',1,'>=',0,'<=',90);
						md = checkfield(md,'fieldname','smb.cotValue','timeseries',1,'NaN',1,'Inf',1,'>=',0);
						md = checkfield(md,'fieldname','smb.ccsnowValue','timeseries',1,'NaN',1,'Inf',1,'>=',0);
						md = checkfield(md,'fieldname','smb.cciceValue','timeseries',1,'NaN',1,'Inf',1,'>=',0);
					end
				case 3
					md = checkfield(md,'fieldname','smb.cldFrac','NaN',1,'Inf',1,'>=',0,'<=',1);
				case 4
					md = checkfield(md,'fieldname','smb.t0wet','NaN',1,'Inf',1,'>=',15,'<=',21.9);
					md = checkfield(md,'fieldname','smb.t0dry','NaN',1,'Inf',1,'>=',30,'<=',30);
					md = checkfield(md,'fieldname','smb.K','NaN',1,'Inf',1,'>=',7,'<=',7);
			end

			%check zTop is < local thickness:
			he=sum(md.geometry.thickness(md.mesh.elements),2)/size(md.mesh.elements,2);
			if any(he<self.zTop),
				error('SMBgemb consistency check error: zTop should be smaller than local ice thickness');
			end
			md = checkfield(md,'fieldname','smb.steps_per_step','>=',1,'numel',[1]);
			md = checkfield(md,'fieldname','smb.averaging','numel',[1],'values',[0 1 2]);
			md = checkfield(md,'fieldname','smb.requested_outputs','stringrow',1);

		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.smb.model','data',8,'format','Integer');

			WriteData(fid,prefix,'object',self,'class','smb','fieldname','isgraingrowth','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','isalbedo','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','isshortwave','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','isthermal','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','isaccumulation','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','ismelt','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','isdensification','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','isturbulentflux','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','isconstrainsurfaceT','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','isdeltaLWup','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','ismappedforcing','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','isprecipforcingremapped','format','Boolean');

			if self.iscompressedforcing
				writetype='CompressedMat';
			else
				writetype='DoubleMat';
			end

			WriteData(fid,prefix,'object',self,'class','smb','fieldname','Ta','format',writetype,'mattype',2,'timeserieslength',size(self.Ta,1),'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','V','format',writetype,'mattype',2,'timeserieslength',size(self.Ta,1),'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','dswrf','format',writetype,'mattype',2,'timeserieslength',size(self.Ta,1),'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','dswdiffrf','format',writetype,'mattype',2,'timeserieslength',size(self.Ta,1),'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','dlwrf','format',writetype,'mattype',2,'timeserieslength',size(self.Ta,1),'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','P','format',writetype,'mattype',2,'timeserieslength',size(self.Ta,1),'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','eAir','format',writetype,'mattype',2,'timeserieslength',size(self.Ta,1),'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','pAir','format',writetype,'mattype',2,'timeserieslength',size(self.Ta,1),'yts',md.constants.yts);

			WriteData(fid,prefix,'object',self,'class','smb','fieldname','Tmean','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','C','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','Vmean','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','Tz','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','Vz','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','zTop','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','dzTop','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','dzMin','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','zY','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','zMax','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','zMin','format','DoubleMat','mattype',2);

			WriteData(fid,prefix,'object',self,'class','smb','fieldname','aIdx','format','Integer');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','eIdx','format','Integer');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','tcIdx','format','Integer');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','swIdx','format','Integer');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','denIdx','format','Integer');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','dsnowIdx','format','Integer');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','InitDensityScaling','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','ThermoDeltaTScaling','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','outputFreq','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','aSnow','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','aIce','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','cldFrac','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','t0wet','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','t0dry','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','K','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','adThresh','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','teThresh','format','Double');

			WriteData(fid,prefix,'object',self,'class','smb','fieldname','aValue','format','DoubleMat','mattype',2,'timeserieslength',md.mesh.numberofelements+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','teValue','format','DoubleMat','mattype',2,'timeserieslength',md.mesh.numberofelements+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','dulwrfValue','format','DoubleMat','mattype',2,'timeserieslength',md.mesh.numberofelements+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','szaValue','format','DoubleMat','mattype',2,'timeserieslength',md.mesh.numberofelements+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','cotValue','format','DoubleMat','mattype',2,'timeserieslength',md.mesh.numberofelements+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','ccsnowValue','format','DoubleMat','mattype',2,'timeserieslength',md.mesh.numberofelements+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','cciceValue','format','DoubleMat','mattype',2,'timeserieslength',md.mesh.numberofelements+1,'yts',md.constants.yts);

			%snow properties init
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','Dzini','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','Dini','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','Reini','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','Gdnini','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','Gspini','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','ECini','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','Wini','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','Aini','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','Adiffini','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','Tini','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','Sizeini','format','IntMat','mattype',2);
			WriteData(fid, prefix, 'object', self, 'fieldname', 'steps_per_step', 'format', 'Integer');
			WriteData(fid, prefix, 'object', self, 'fieldname', 'averaging', 'format', 'Integer');

			if (self.ismappedforcing)
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','mappedforcingpoint','format','IntMat','mattype',2);
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','mappedforcingelevation','format','DoubleMat','mattype',3);
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','lapseTaValue','format','DoubleMat','mattype',3);
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','lapsedlwrfValue','format','DoubleMat','mattype',3);
			end

			%figure out dt from forcings:
			if (any(md.smb.P(end,:) - md.smb.Ta(end,:) ~= 0) | ...
					any(md.smb.V(end,:) - md.smb.Ta(end,:) ~= 0) | ...
					any(md.smb.dswrf(end,:) - md.smb.Ta(end,:) ~= 0) | ...
					any(md.smb.dlwrf(end,:) - md.smb.Ta(end,:) ~= 0) | ...
					any(md.smb.eAir(end,:) - md.smb.Ta(end,:) ~= 0) | ...
					any(md.smb.pAir(end,:) - md.smb.Ta(end,:) ~= 0))
				error('All GEMB forcings (Ta, P, V, dswrf, dlwrf, eAir, pAir) must have the same time steps in the final row!');
			end
			if size(md.smb.teValue,2)>1 & any(md.smb.teValue(end,:) - md.smb.Ta(end,:) ~= 0)
				error('If GEMB forcing teValue is transient, it must have the same time steps as input Ta in the final row!');
			end
			if size(md.smb.dswdiffrf,2)>1 & any(md.smb.dswdiffrf(end,:) - md.smb.Ta(end,:) ~= 0)
				error('If GEMB forcing dswdiffrf is transient, it must have the same time steps as input Ta in the final row!');
			end
			if size(md.smb.dulwrfValue,2)>1 & any(md.smb.dulwrfValue(end,:) - md.smb.Ta(end,:) ~= 0)
				error('If GEMB forcing dulwrfValue is transient, it must have the same time steps as input Ta in the final row!');
			end
			if size(md.smb.aValue,2)>1 & any(md.smb.aValue(end,:) - md.smb.Ta(end,:) ~= 0)
				error('If GEMB forcing aValue is transient, it must have the same time steps as input Ta in the final row!');
			end
			if size(md.smb.szaValue,2)>1 & any(md.smb.szaValue(end,:) - md.smb.Ta(end,:) ~= 0)
				error('If GEMB forcing szaValue is transient, it must have the same time steps as input Ta in the final row!');
			end
			if size(md.smb.cotValue,2)>1 & any(md.smb.cotValue(end,:) - md.smb.Ta(end,:) ~= 0)
				error('If GEMB forcing cotValue is transient, it must have the same time steps as input Ta in the final row!');
			end
			if size(md.smb.ccsnowValue,2)>1 & any(md.smb.ccsnowValue(end,:) - md.smb.Ta(end,:) ~= 0)
				error('If GEMB forcing ccsnowValue is transient, it must have the same time steps as input Ta in the final row!');
			end
			if size(md.smb.cciceValue,2)>1 & any(md.smb.cciceValue(end,:) - md.smb.Ta(end,:) ~= 0)
				error('If GEMB forcing cciceValue is transient, it must have the same time steps as input Ta in the final row!');
			end
			time=self.Ta(end,:); %assume all forcings are on the same time step
			dtime=diff(time,1);
			dt=min(dtime);

			WriteData(fid,prefix,'data',dt,'name','md.smb.dt','format','Double','scale',yts);

			% Check if smb_dt goes evenly into transient core time step
			if (mod(md.timestepping.time_step,dt) >= 1e-10)
				error('smb_dt/dt = %f. The number of SMB time steps in one transient core time step has to be an an integer',md.timestepping.time_step/dt);
			end

			% Make sure that adaptive time step is off
			if isa(md.timestepping,'timesteppingadaptive')
				error('GEMB cannot be run with adaptive timestepping.  Check class type of md.timestepping');
			end

			%process requested outputs
			outputs = self.requested_outputs;
			pos  = find(ismember(outputs,'default'));
			if ~isempty(pos),
				outputs(pos) = [];                         %remove 'default' from outputs
				outputs      = [outputs defaultoutputs(self,md)]; %add defaults
			end
			WriteData(fid,prefix,'data',outputs,'name','md.smb.requested_outputs','format','StringArray');
		end % }}}
	end
end
