%MATICE class definition
%
%   Usage:
%      matice=matice();

classdef matice
	properties (SetAccess=public) 
		rho_ice                         = 0.;
		rho_water                       = 0.;
		rho_freshwater                  = 0.;
		mu_water                        = 0.;
		heatcapacity                    = 0.;
		latentheat                      = 0.;
		thermalconductivity             = 0.;
		temperateiceconductivity        = 0.;
		effectiveconductivity_averaging = 0.;
		meltingpoint                    = 0.;
		beta                            = 0.;
		mixed_layer_capacity            = 0.;
		thermal_exchange_velocity       = 0.;
		rheology_B                      = NaN;
		rheology_n                      = NaN;
		rheology_law                    = '';

		%slc
		earth_density                   = 0;

	end
	methods
		function self = extrude(self,md) % {{{
			self.rheology_B=project3d(md,'vector',self.rheology_B,'type','node');
			self.rheology_n=project3d(md,'vector',self.rheology_n,'type','element');
		end % }}}
		function self = matice(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					inputstruct=varargin{1};
					list1 = properties('matice');
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
		function self = setdefaultparameters(self) % {{{

			%ice density (kg/m^3)
			self.rho_ice=917.;

			%ocean water density (kg/m^3)
			self.rho_water=1023.;

			%fresh water density (kg/m^3)
			self.rho_freshwater=1000.;

			%water viscosity (N.s/m^2)
			self.mu_water=0.001787;

			%ice heat capacity cp (J/kg/K)
			self.heatcapacity=2093.;

			%ice latent heat of fusion L (J/kg)
			self.latentheat=3.34*10^5;

			%ice thermal conductivity (W/m/K)
			self.thermalconductivity=2.4;
			
			%wet ice thermal conductivity (W/m/K)
			self.temperateiceconductivity=.24;

			%computation of effective conductivity
			self.effectiveconductivity_averaging=1;

			%the melting point of ice at 1 atmosphere of pressure in K
			self.meltingpoint=273.15;

			%rate of change of melting point with pressure (K/Pa)
			self.beta=9.8*10^-8;

			%mixed layer (ice-water interface) heat capacity (J/kg/K)
			self.mixed_layer_capacity=3974.;

			%thermal exchange velocity (ice-water interface) (m/s)
			self.thermal_exchange_velocity=1.00*10^-4;

			%Rheology law: what is the temperature dependence of B with T
			%available: none, paterson and arrhenius
			self.rheology_law='Paterson';

			%Rheology for ice: 
			self.rheology_B=2.1*10^8;
			self.rheology_n=3;

			%SLR
			self.earth_density= 5512;  % average density of the Earth, (kg/m^3)

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if strcmpi(solution,'TransientSolution') & md.transient.isslc,
				md = checkfield(md,'fieldname','materials.earth_density','>',0,'numel',1);
			else
				md = checkfield(md,'fieldname','materials.rho_ice','>',0);
				md = checkfield(md,'fieldname','materials.rho_water','>',0);
				md = checkfield(md,'fieldname','materials.rho_freshwater','>',0);
				md = checkfield(md,'fieldname','materials.mu_water','>',0);
				md = checkfield(md,'fieldname','materials.rheology_B','>',0,'size','universal','NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','materials.rheology_n','>',0,'size','universal','NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','materials.rheology_law','values',{'None' 'BuddJacka' 'Cuffey' 'CuffeyTemperate' 'Paterson' 'Arrhenius' 'LliboutryDuval' 'NyeCO2' 'NyeH2O'});
				md = checkfield(md,'fieldname','materials.effectiveconductivity_averaging','numel',[1],'values',[0 1 2]);
			end

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   Materials:'));

			fielddisplay(self,'rho_ice','ice density [kg/m^3]');
			fielddisplay(self,'rho_water','ocean water density [kg/m^3]');
			fielddisplay(self,'rho_freshwater','fresh water density [kg/m^3]');
			fielddisplay(self,'mu_water','water viscosity [N s/m^2]');
			fielddisplay(self,'heatcapacity','heat capacity [J/kg/K]');
			fielddisplay(self,'thermalconductivity',['ice thermal conductivity [W/m/K]']);
			fielddisplay(self,'temperateiceconductivity','temperate ice thermal conductivity [W/m/K]');
			fielddisplay(self,'effectiveconductivity_averaging','computation of effective conductivity: (0) arithmetic mean, (1) harmonic mean, (2) geometric mean (default)');
			fielddisplay(self,'meltingpoint','melting point of ice at 1atm in K');
			fielddisplay(self,'latentheat','latent heat of fusion [J/kg]');
			fielddisplay(self,'beta','rate of change of melting point with pressure [K/Pa]');
			fielddisplay(self,'mixed_layer_capacity','mixed layer capacity [W/kg/K]');
			fielddisplay(self,'thermal_exchange_velocity','thermal exchange velocity [m/s]');
			fielddisplay(self,'rheology_B','flow law parameter [Pa s^(1/n)]');
			fielddisplay(self,'rheology_n','Glen''s flow law exponent');
			fielddisplay(self,'rheology_law',['law for the temperature dependance of the rheology: ''None'', ''BuddJacka'', Cuffey'', ''CuffeyTemperate'', ''Paterson'', ''Arrhenius'', ''LliboutryDuval'', ''NyeH2O'', or ''NyeCO2''']);
			fielddisplay(self,'earth_density','Mantle density [kg/m^-3]');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'name','md.materials.type','data',3,'format','Integer');
			WriteData(fid,prefix,'object',self,'class','materials','fieldname','rho_ice','format','Double');
			WriteData(fid,prefix,'object',self,'class','materials','fieldname','rho_water','format','Double');
			WriteData(fid,prefix,'object',self,'class','materials','fieldname','rho_freshwater','format','Double');
			WriteData(fid,prefix,'object',self,'class','materials','fieldname','mu_water','format','Double');
			WriteData(fid,prefix,'object',self,'class','materials','fieldname','heatcapacity','format','Double');
			WriteData(fid,prefix,'object',self,'class','materials','fieldname','latentheat','format','Double');
			WriteData(fid,prefix,'object',self,'class','materials','fieldname','thermalconductivity','format','Double');
			WriteData(fid,prefix,'object',self,'class','materials','fieldname','temperateiceconductivity','format','Double');
			WriteData(fid,prefix,'object',self,'class','materials','fieldname','effectiveconductivity_averaging','format','Integer');
			WriteData(fid,prefix,'object',self,'class','materials','fieldname','meltingpoint','format','Double');
			WriteData(fid,prefix,'object',self,'class','materials','fieldname','beta','format','Double');
			WriteData(fid,prefix,'object',self,'class','materials','fieldname','mixed_layer_capacity','format','Double');
			WriteData(fid,prefix,'object',self,'class','materials','fieldname','thermal_exchange_velocity','format','Double');
			if(size(self.rheology_B,1)==md.mesh.numberofvertices | size(self.rheology_B,1)==md.mesh.numberofvertices+1 | (size(self.rheology_B,1)==md.mesh.numberofelements && size(self.rheology_B,2)>1))
				mattype=1; tsl = md.mesh.numberofvertices;
			else
				mattype=2; tsl = md.mesh.numberofelements;
			end
			WriteData(fid,prefix,'object',self,'class','materials','fieldname','rheology_B','format','DoubleMat','mattype',mattype,'timeserieslength',tsl+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','materials','fieldname','rheology_n','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'data',self.rheology_law,'name','md.materials.rheology_law','format','String');
			WriteData(fid,prefix,'object',self,'class','materials','fieldname','earth_density','format','Double');
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejsdouble(fid,[modelname '.materials.rho_ice'],self.rho_ice);
			writejsdouble(fid,[modelname '.materials.rho_water'],self.rho_water);
			writejsdouble(fid,[modelname '.materials.rho_freshwater'],self.rho_freshwater);
			writejsdouble(fid,[modelname '.materials.mu_water'],self.mu_water);
			writejsdouble(fid,[modelname '.materials.heatcapacity'],self.heatcapacity);
			writejsdouble(fid,[modelname '.materials.latentheat'],self.latentheat);
			writejsdouble(fid,[modelname '.materials.thermalconductivity'],self.thermalconductivity);
			writejsdouble(fid,[modelname '.materials.temperateiceconductivity'],self.temperateiceconductivity);
			writejsdouble(fid,[modelname '.materials.effectiveconductivity_averaging'],self.effectiveconductivity_averaging);
			writejsdouble(fid,[modelname '.materials.meltingpoint'],self.meltingpoint);
			writejsdouble(fid,[modelname '.materials.beta'],self.beta);
			writejsdouble(fid,[modelname '.materials.mixed_layer_capacity'],self.mixed_layer_capacity);
			writejsdouble(fid,[modelname '.materials.thermal_exchange_velocity'],self.thermal_exchange_velocity);
			writejsdouble(fid,[modelname '.materials.mixed_layer_capacity'],self.mixed_layer_capacity);
			writejs1Darray(fid,[modelname '.materials.rheology_B'],self.rheology_B);
			writejs1Darray(fid,[modelname '.materials.rheology_n'],self.rheology_n);
			writejsstring(fid,[modelname '.materials.rheology_law'],self.rheology_law);
			writejsdouble(fid,[modelname '.materials.earth_density'],self.earth_density);

		end % }}}
	end
end
