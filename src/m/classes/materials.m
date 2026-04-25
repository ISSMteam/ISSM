%MATERIALS class definition
%
%   Usage:
%      materials=materials();

classdef materials < dynamicprops
	properties (SetAccess=public)
		nature={};
		%all properties are dynamic.
	end
	methods
		function self = materials(varargin) % {{{
			if nargin==0
				self.nature={'ice'};
			else
				self.nature=varargin;
			end

			%check this is acceptable:
			for i=1:length(self.nature),
				if ~(strcmpi(self.nature{i},'litho') | strcmpi(self.nature{i},'ice') | strcmpi(self.nature{i},'hydro')),
					error('materials constructor error message: nature of the material not supported yet! (''ice'' or ''litho'' or ''hydro'')');
				end
			end

			%start filling in the dynamic fields:
			for i=1:length(self.nature),
				nat=self.nature{i};
				switch nat
				case 'ice'
					self.addprop('rho_ice');
					self.addprop('rho_water');
					self.addprop('rho_freshwater');
					self.addprop('mu_water');
					self.addprop('heatcapacity');
					self.addprop('latentheat');
					self.addprop('thermalconductivity');
					self.addprop('temperateiceconductivity');
					self.addprop('effectiveconductivity_averaging');
					self.addprop('meltingpoint');
					self.addprop('beta');
					self.addprop('mixed_layer_capacity');
					self.addprop('thermal_exchange_velocity');
					self.addprop('rheology_B');
					self.addprop('rheology_n');
					self.addprop('rheology_law');
				case 'litho'
					self.addprop('numlayers');
					self.addprop('radius');
					self.addprop('viscosity');
					self.addprop('lame_lambda');
					self.addprop('lame_mu');
					self.addprop('burgers_viscosity');
					self.addprop('burgers_mu');
					self.addprop('ebm_alpha');
					self.addprop('ebm_delta');
					self.addprop('ebm_taul');
					self.addprop('ebm_tauh');
					self.addprop('rheologymodel');
					self.addprop('density');
					self.addprop('issolid');
				case 'hydro'
					self.addprop('rho_ice');
					self.addprop('rho_water');
					self.addprop('rho_freshwater');
				otherwise
					error('materials constructor error message: nature of the material not supported yet! (''ice'' or ''litho'' or ''hydro'')');
				end
			end
			self.addprop('earth_density');

			%set default parameters:
			self.setdefaultparameters();

		end % }}}
		function self = setdefaultparameters(self) % {{{

			for i=1:length(self.nature),
				nat=self.nature{i};
				switch nat
				case 'ice'
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
					self.latentheat=3.34*1e5;

					%ice thermal conductivity (W/m/K)
					self.thermalconductivity=2.4;

					%wet ice thermal conductivity (W/m/K)
					self.temperateiceconductivity=.24;

					%computation of effective conductivity
					self.effectiveconductivity_averaging=1;

					%the melting point of ice at 1 atmosphere of pressure in K
					self.meltingpoint=273.15;

					%rate of change of melting point with pressure (K/Pa)
					self.beta=9.8*1e-8;

					%mixed layer (ice-water interface) heat capacity (J/kg/K)
					self.mixed_layer_capacity=3974.;

					%thermal exchange velocity (ice-water interface) (m/s)
					self.thermal_exchange_velocity=1.00*1e-4;

					%Rheology law: what is the temperature dependence of B with T
					%available: none, paterson and arrhenius
					self.rheology_law='Paterson';

					%Rheology fields default: 
					self.rheology_B   = 1 * 1e8;
					self.rheology_n   = 3;

				case 'litho'
					%we default to a configuration that enables running GIA solutions using giacaron and/or giaivins.
					self.numlayers=2;

					%center of the earth (approximation, must not be 0), then the lab (lithosphere/asthenosphere boundary) then the surface
					%(with 1d3 to avoid numerical singularities)
					self.radius=[1e3;6278*1e3;6378*1e3];

					self.viscosity=[1e21;1e40]; %mantle and lithosphere viscosity (respectively) [Pa.s]
					self.lame_mu=[1.45*1e11;6.7*1e10];  % (Pa) %lithosphere and mantle shear modulus (respectively) [Pa]
					self.lame_lambda=self.lame_mu;  % (Pa) %mantle and lithosphere lamba parameter (respectively) [Pa]
					self.burgers_viscosity=[NaN;NaN];
					self.burgers_mu=[NaN;NaN];

					self.ebm_alpha=[NaN;NaN];
					self.ebm_delta=[NaN;NaN];
					self.ebm_taul=[NaN;NaN];
					self.ebm_tauh=[NaN;NaN];
					self.rheologymodel=[0;0];
					self.density=[5.51*1e3;5.50*1e3];  % (Pa) %mantle and lithosphere density [kg/m^3]
					self.issolid=[1;1]; % is layer solid or liquid.

				case 'hydro'
					%ice density (kg/m^3)
					self.rho_ice=917.;

					%ocean water density (kg/m^3)
					self.rho_water=1023.;
					
					%fresh water density (kg/m^3)
					self.rho_freshwater=1000.;

				otherwise
					error('materials constructor error message: nature of the material not supported yet! (''ice'' or ''litho'' or ''hydro'')');
				end

				% average density of the Earth (kg/m^3)
				self.earth_density=5512;

			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   Materials:'));

			for i=1:length(self.nature),
				nat=self.nature{i};
				switch nat
				case 'ice'
					disp(sprintf('\n      Ice:'));
					fielddisplay(self,'rho_ice','ice density [kg/m^3]');
					fielddisplay(self,'rho_water','ocean water density [kg/m^3]');
					fielddisplay(self,'rho_freshwater','fresh water density [kg/m^3]');
					fielddisplay(self,'mu_water','water viscosity [N s/m^2]');
					fielddisplay(self,'heatcapacity','heat capacity [J/kg/K]');
					fielddisplay(self,'thermalconductivity',['ice thermal conductivity [W/m/K]']);
					fielddisplay(self,'temperateiceconductivity','temperate ice thermal conductivity [W/m/K]');
					fielddisplay(self,'meltingpoint','melting point of ice at 1atm in K');
					fielddisplay(self,'latentheat','latent heat of fusion [J/kg]');
					fielddisplay(self,'beta','rate of change of melting point with pressure [K/Pa]');
					fielddisplay(self,'mixed_layer_capacity','mixed layer capacity [W/kg/K]');
					fielddisplay(self,'thermal_exchange_velocity','thermal exchange velocity [m/s]');
					fielddisplay(self,'rheology_B','flow law parameter [Pa s^(1/n)]');
					fielddisplay(self,'rheology_n','Glen''s flow law exponent');
					fielddisplay(self,'rheology_law',['law for the temperature dependance of the rheology: ''None'', ''BuddJacka'', Cuffey'', ''CuffeyTemperate'', ''Paterson'', ''Arrhenius'', ''LliboutryDuval'', ''NyeCO2'', or ''NyeH2O''']);
				case 'litho'
					disp(sprintf('\n      Litho:'));
					fielddisplay(self,'numlayers','number of layers (default: 2)');
					fielddisplay(self,'radius','array describing the radius for each interface (numlayers+1) [m]');
					fielddisplay(self,'viscosity','array describing each layer''s viscosity (numlayers) [Pa.s]');
					fielddisplay(self,'lame_lambda','array describing the lame lambda parameter (numlayers) [Pa]');
					fielddisplay(self,'lame_mu','array describing the shear modulus for each layers (numlayers) [Pa]');
					fielddisplay(self,'burgers_viscosity','array describing each layer''s transient viscosity, only for Burgers rheologies  (numlayers) [Pa.s]');
					fielddisplay(self,'burgers_mu','array describing each layer''s transient shear modulus, only for Burgers rheologies  (numlayers) [Pa]');

					fielddisplay(self,'ebm_alpha','array describing each layer''s exponent parameter controlling the shape of shear modulus curve between taul and tauh, only for EBM rheology (numlayers)');
					fielddisplay(self,'ebm_delta','array describing each layer''s amplitude of the transient relaxation (ratio between elastic rigity to pre-maxwell relaxation rigity), only for EBM rheology (numlayers)');
					fielddisplay(self,'ebm_taul','array describing each layer''s starting period for transient relaxation, only for EBM rheology  (numlayers) [s]');
					fielddisplay(self,'ebm_tauh','array describing each layer''s array describing each layer''s end period for transient relaxation, only for Burgers rheology (numlayers) [s]');


					fielddisplay(self,'rheologymodel','array describing whether we adopt a Maxwell (0), Burgers (1) or EBM (2) rheology (default: 0)');
					fielddisplay(self,'density','array describing each layer''s density (numlayers) [kg/m^3]');
					fielddisplay(self,'issolid','array describing whether the layer is solid or liquid (default 1) (numlayers)');
				case 'hydro'
					disp(sprintf('\n      Hydro:'));
					fielddisplay(self,'rho_ice','ice density [kg/m^3]');
					fielddisplay(self,'rho_water','ocean water density [kg/m^3]');
					fielddisplay(self,'earth_density','mantle density [kg/m^3]');
					fielddisplay(self,'rho_freshwater','fresh water density [kg/m^3]');

				otherwise
					error('materials constructor error message: nature of the material not supported yet! (''ice'' or ''litho'' or ''hydro'')');
				end
			end
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			for i=1:length(self.nature),
				nat=self.nature{i};
				switch nat
				case 'ice'
					md = checkfield(md,'fieldname','materials.rho_ice','>',0);
					md = checkfield(md,'fieldname','materials.rho_water','>',0);
					md = checkfield(md,'fieldname','materials.rho_freshwater','>',0);
					md = checkfield(md,'fieldname','materials.mu_water','>',0);
					md = checkfield(md,'fieldname','materials.rheology_B','>',0,'timeseries',1,'NaN',1,'Inf',1);
					md = checkfield(md,'fieldname','materials.rheology_n','>',0,'size',[md.mesh.numberofelements 1]);
					md = checkfield(md,'fieldname','materials.rheology_law','values',{'None' 'BuddJacka' 'Cuffey' 'CuffeyTemperate' 'Paterson' 'Arrhenius' 'LliboutryDuval' 'NyeCO2' 'NyeH2O'});
				case 'litho'
					if ~ismember('LoveAnalysis',analyses), return; end
					md = checkfield(md,'fieldname','materials.numlayers','NaN',1,'Inf',1,'>',0,'numel',1);
					md = checkfield(md,'fieldname','materials.radius','NaN',1,'Inf',1,'size',[md.materials.numlayers+1 1],'>',0);
					md = checkfield(md,'fieldname','materials.lame_mu','NaN',1,'Inf',1,'size',[md.materials.numlayers 1],'>=',0);
					md = checkfield(md,'fieldname','materials.lame_lambda','NaN',1,'Inf',1,'size',[md.materials.numlayers 1],'>=',0);
					md = checkfield(md,'fieldname','materials.issolid','NaN',1,'Inf',1,'size',[md.materials.numlayers 1],'>=',0,'<',2);
					md = checkfield(md,'fieldname','materials.density','NaN',1,'Inf',1,'size',[md.materials.numlayers 1],'>',0);
					md = checkfield(md,'fieldname','materials.viscosity','NaN',1,'Inf',1,'size',[md.materials.numlayers 1],'>=',0);
					md = checkfield(md,'fieldname','materials.rheologymodel','NaN',1,'Inf',1,'size',[md.materials.numlayers 1],'>=',0,'<=',2);
					if any(self.rheologymodel==1)
						md = checkfield(md,'fieldname','materials.burgers_viscosity','Inf',1,'size',[md.materials.numlayers 1],'>=',0);
						md = checkfield(md,'fieldname','materials.burgers_mu','Inf',1,'size',[md.materials.numlayers 1],'>=',0);
					end
					if any(self.rheologymodel==2)
						md = checkfield(md,'fieldname','materials.ebm_alpha','Inf',1,'size',[md.materials.numlayers 1],'>=',0);
						md = checkfield(md,'fieldname','materials.ebm_delta','Inf',1,'size',[md.materials.numlayers 1],'>=',0);
						md = checkfield(md,'fieldname','materials.ebm_taul','Inf',1,'size',[md.materials.numlayers 1],'>=',0);
						md = checkfield(md,'fieldname','materials.ebm_tauh','Inf',1,'size',[md.materials.numlayers 1],'>=',0);
					end
                    if any(diff(md.materials.radius)<=0)
							error('materials checkconsistency error message: radius should be monotonously increasing');                        
                    end
					for i=1:md.materials.numlayers,
						if md.materials.rheologymodel(i)==1 & (isnan(md.materials.burgers_viscosity(i) | isnan(md.materials.burgers_mu(i)))),
							error('materials checkconsistency error message: Litho burgers_viscosity or burgers_mu has NaN values, inconsistent with rheologymodel choice');
						end
						if md.materials.rheologymodel(i)==2 & (isnan(md.materials.ebm_alpha(i)) | isnan(md.materials.ebm_delta(i)) | isnan(md.materials.ebm_taul(i)) | isnan(md.materials.ebm_tauh(i))),
							error('materials checkconsistency error message: Litho ebm_alpha, ebm_delta, ebm_taul or ebm_tauh has NaN values, inconsistent with rheologymodel choice');
						end
					end
					if md.materials.issolid(1)==0 | md.materials.lame_mu(1)==0
						error('First layer must be solid (issolid(1) > 0 AND lame_mu(1) > 0). Add a weak inner core if necessary.');
					end
					ind=find(md.materials.issolid==0);
					%if sum(ismember(diff(ind),1)>=1) %if there are at least two consecutive indices that contain issolid = 0
					%		error(['Fluid layers detected at layers #', num2str(ind'), ', but having 2 or more adjacent fluid layers is not supported yet. Consider merging them.'])
					%end

				case 'hydro'
					md = checkfield(md,'fieldname','materials.rho_ice','>',0);
					md = checkfield(md,'fieldname','materials.rho_water','>',0);
					md = checkfield(md,'fieldname','materials.earth_density','>',0,'numel',1);
					md = checkfield(md,'fieldname','materials.rho_freshwater','>',0);

				otherwise
					error('materials constructor error message: nature of the material not supported yet! (''ice'' or ''litho'' or ''hydro'')');
				end
			end

		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			%1: MatdamageiceEnum 2: MatestarEnum 3: MaticeEnum 4: MatenhancediceEnum 5: MaterialsEnum
			WriteData(fid,prefix,'name','md.materials.nature','data',naturetointeger(self.nature),'format','IntMat','mattype',3);
			WriteData(fid,prefix,'name','md.materials.type','data',5,'format','Integer'); %DANGER: this can evolve if you have classes.
			for i=1:length(self.nature),
				nat=self.nature{i};
				switch nat
				case 'ice'
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
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','rheology_B','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','rheology_n','format','DoubleMat','mattype',2);
					WriteData(fid,prefix,'data',self.rheology_law,'name','md.materials.rheology_law','format','String');
				case 'litho'
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','numlayers','format','Integer');
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','radius','format','DoubleMat','mattype',3);
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','lame_mu','format','DoubleMat','mattype',3);
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','lame_lambda','format','DoubleMat','mattype',3);
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','issolid','format','DoubleMat','mattype',3);
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','density','format','DoubleMat','mattype',3);
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','viscosity','format','DoubleMat','mattype',3);
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','rheologymodel','format','DoubleMat','mattype',3);
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','burgers_viscosity','format','DoubleMat','mattype',3);
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','burgers_mu','format','DoubleMat','mattype',3);
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','ebm_alpha','format','DoubleMat','mattype',3);
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','ebm_delta','format','DoubleMat','mattype',3);
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','ebm_taul','format','DoubleMat','mattype',3);
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','ebm_tauh','format','DoubleMat','mattype',3);
					%compute earth density compatible with our layer density distribution: 
					earth_density=0;
					for i=1:self.numlayers,
						earth_density=earth_density + (self.radius(i+1)^3-self.radius(i)^3)*self.density(i);
					end
					earth_density=earth_density/self.radius(self.numlayers+1)^3;
					self.earth_density=earth_density;
				case 'hydro'
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','rho_ice','format','Double');
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','rho_water','format','Double');
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','rho_freshwater','format','Double');
				otherwise
					error('materials constructor error message: nature of the material not supported yet! (''ice'' or ''litho'' or ''hydro'')');
				end
			end
			WriteData(fid,prefix,'data',self.earth_density,'name','md.materials.earth_density','format','Double');
		end % }}}
		function self = extrude(self,md) % {{{
			for i=1:length(self.nature),
				nat=self.nature{i};
				switch nat
				case 'ice'
					self.rheology_B=project3d(md,'vector',self.rheology_B,'type','node');
					self.rheology_n=project3d(md,'vector',self.rheology_n,'type','element');
				end
			end
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{

			for i=1:length(self.nature),
				nat=self.nature{i};
				switch nat
				case 'ice'
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
				case 'litho'
					writejsdouble(fid,[modelname '.materials.numlayers'],self.numlayers);
					writejsdouble(fid,[modelname '.materials.radius'],self.radius);
					writejsdouble(fid,[modelname '.materials.lame_mu'],self.lame_mu);
					writejsdouble(fid,[modelname '.materials.lame_lambda'],self.lame_lambda);
					writejsdouble(fid,[modelname '.materials.issolid'],self.issolid);
					writejsdouble(fid,[modelname '.materials.density'],self.density);
					writejsdouble(fid,[modelname '.materials.viscosity'],self.viscosity);
					writejsdouble(fid,[modelname '.materials.rheologymodel'],self.rheologymodel);
					writejsdouble(fid,[modelname '.materials.burgers_viscosity'],self.burgers_viscosity);
					writejsdouble(fid,[modelname '.materials.burgers_mu'],self.burgers_mu);
					writejsdouble(fid,[modelname '.materials.ebm_alpha'],self.ebm_alpha);
					writejsdouble(fid,[modelname '.materials.ebm_delta'],self.ebm_delta);
					writejsdouble(fid,[modelname '.materials.ebm_taul'],self.ebm_taul);
					writejsdouble(fid,[modelname '.materials.ebm_tauh'],self.ebm_tauh);

				case 'hydro'
					writejsdouble(fid,[modelname '.materials.rho_ice'],self.rho_ice);
					writejsdouble(fid,[modelname '.materials.rho_water'],self.rho_water);
					writejsdouble(fid,[modelname '.materials.earth_density'],self.earth_density);
					writejsdouble(fid,[modelname '.materials.rho_freshwater'],self.rho_freshwater);

				otherwise
					error('materials constructor error message: nature of the material not supported yet! (''ice'' or ''litho'' or ''hydro'')');
				end
			end



		end % }}}
		function setlitho2prem(self) % {{{
			%based on materials.radius and materials.issolid, replaces materials.density, materials.lame_mu, materials.lame_lambda
			% with volumetric averages of a polynomial version of the PREM model
			if ~strcmpi(self.nature, 'litho')
			error('materials setlitho2prem error message: materials.nature must be ''litho''')
			end
		   
		   
		   ra=self.radius(end);
		   r(1) = 0.0;     r(2) = 1221.5; r(3) = 3480.0; r(4) = 3630.;
		   r(5) = 5600.0;  r(6) = 5701.0;  r(7) = 5771.0; r(8) = 5971.;
		   r(9) = 6151.0; r(10) = 6291.0; r(11) = 6346.6;
		   r(12) = 6356.0; r(13) = 6368.0; r(14) = 6371.;
		  
		   if (r(14)*1e3 ~= ra) 
			error(['materials setlitho2prem error message: Earth surface radius in PREM r=' num2str(r(end)*1e3) ' does not match materials.radius(end)=' num2str(self.radius(end))]); 
		   end

		  d = zeros(12,4); %polynomial coef for density
		  d(1,1) = 13.0885; 		         	  d(1,3) = -8.8381 ;
		  d(2,1) = 12.5815; d(2,2) = -1.2638; d(2,3) = -3.6426; d(2,4) = -5.5281;
		  d(3,1) = 7.9565 ; d(3,2) = -6.4761;   d(3,3) = 5.5283;  d(3,4) = -3.0807;
		  d(4,1) = 7.9565 ; d(4,2) = -6.4761;   d(4,3) = 5.5283;  d(4,4) = -3.0807;
		  d(5,1) = 7.9565 ; d(5,2) = -6.4761;   d(5,3) = 5.5283;  d(5,4) = -3.0807;
		  d(6,1) = 5.3197 ; d(6,2) = -1.4836;
		  d(7,1) = 11.2494; d(7,2) = -8.0298;
		  d(8,1) = 7.1089 ; d(8,2) = -3.8045;
		  d(9,1) = 2.6910 ; d(9,2) = 0.6924;
		  d(10,1) = 2.6910; d(10,2) = 0.6924;
		  d(11,1) = 2.9  ;
		  d(12,1) = 2.6  ;

		% ocean
		  if (~self.issolid(end)) 
		  	d(13,1) = 1.02 ;

		% continental
		  else
		  d(13,1) = d(12,1);
		  end

		  p = zeros(13,4); %polynomial coef for Vp, pressure wave velocity
		  p(1,1) = 11.2622 ; p(1,3) = -6.3640;
		  p(2,1) = 11.0487 ; p(2,2) = -4.0362; p(2,3)  = 4.8023; p(2,4) = -13.5732;
		  p(3,1) = 15.3891 ; p(3,2) = -5.3181; p(3,3)  = 5.5242; p(3,4) = -2.5514;
		  p(4,1) = 24.952 ; p(4,2)  = -40.4673; p(4,3) = 51.4832; p(4,4) = -26.6419;
		  p(5,1) = 29.2766 ; p(5,2) = -23.6027; p(5,3) = 5.5242; p(5,4) = -2.5514;
		  p(6,1) = 19.0957 ; p(6,2)  = -9.8672;
		  p(7,1) = 39.7027 ; p(7,2)  = -32.6166;
		  p(8,1) = 20.3926 ; p(8,2)  = -12.2569;
		  p(9,1) = 4.1875 ; p(9,2)  = 3.9382;
		  p(10,1) = 4.1875 ; p(10,2) = 3.9382;
		  p(11,1) = 6.8 ;
		  p(12,1) = 5.8;
		%
		% ocean
		  if (~self.issolid(end)) 
		  p(13,1) = 1.45 ;
		%
		% continental
		  else
		  p(13,1) = p(12,1);
		  end
		%----
		%
		  s = zeros(13,4); %polynomial coef for Vs, shear wave velocity
		%
		  s(1,1) = 3.6678; s(1,3) = -4.4475;

		  s(3,1) = 6.9254; s(3,2) = 1.4672; s(3,3) = -2.0834; s(3,4) = 0.9783;
		  s(4,1) = 11.1671; s(4,2) = -13.7818; s(4,3) = 17.4575; s(4,4) = -9.2777;
		  s(5,1) = 22.3459; s(5,2) = -17.2473; s(5,3) = -2.0834; s(5,4) = 0.9783;
		  s(6,1) = 9.9839; s(6,2) = -4.9324;
		  s(7,1) = 22.3512; s(7,2) = -18.5856 ;
		  s(8,1) = 8.9496; s(8,2) = -4.4597;
		  s(9,1) = 2.1519; s(9,2) = 2.3481;
		  s(10,1) = 2.1519; s(10,2) = 2.3481;
		  s(11,1) = 3.9 ;
		  s(12,1) = 3.2 ;
		%
		% ocean (please don't modify)
		  if (~self.issolid(end))
		%
		% continental
		  else
		  s(13,1) = s(12,1);
		  end
		%
		%
		  r = r*1e3;
		  
		  %- handling the first layer : central sphere
		  rad = self.radius;
		  rad(1) = 0.;
		  
		  for j = 1:self.numlayers
		  
			ro = 0.;
			vp = 0.;
			vs = 0.;

			for i = 1:13
				r1 = 0.;
				r2 = 0.;
				if ((rad(j) > r(i)) & (rad(j) <= r(i+1))) 
					if (rad(j+1) <= r(i+1)) 
						r2 = rad(j+1);
						r1 = rad(j);
					else
						r2 = r(i+1);
						r1 = rad(j);
					end
				elseif (rad(j) <= r(i))
					if ((rad(j+1) > r(i)) & (rad(j+1) <= r(i+1)))
						r2 = rad(j+1);
						r1 = r(i);
					elseif (rad(j+1) > r(i+1))
						r2 = r(i+1);
						r1 = r(i);
					end
				end

				t1 = d(i,1)/3.;
				t2 = d(i,2)/(ra*4.);
				t3 = d(i,3)/((ra^2)*5.);
				t4 = d(i,4)/((ra^3)*6.);
				ro =  ro + t1*(r2^3) + t2*(r2^4) + t3*(r2^5) + t4*(r2^6) - ( t1*(r1^3) + t2*(r1^4) + t3*(r1^5) + t4*(r1^6) );
					  
				t1 = p(i,1)/3.;
				t2 = p(i,2)/(ra*4.);
				t3 = p(i,3)/((ra^2)*5.);
				t4 = p(i,4)/((ra^3)*6.);
				vp =  vp + t1*(r2^3) + t2*(r2^4) + t3*(r2^5) + t4*(r2^6) - ( t1*(r1^3) + t2*(r1^4) + t3*(r1^5) + t4*(r1^6) );
					  
				t1 = s(i,1)/3.;
				t2 = s(i,2)/(ra*4.);
				t3 = s(i,3)/((ra^2)*5.);
				t4 = s(i,4)/((ra^3)*6.);
				vs =  vs + t1*(r2^3) + t2*(r2^4) + t3*(r2^5) + t4*(r2^6) - ( t1*(r1^3) + t2*(r1^4) + t3*(r1^5) + t4*(r1^6) );

			end
			ro = ro*3 / (rad(j+1)^3-rad(j)^3);
			vp = vp*3 /(rad(j+1)^3-rad(j)^3);
			vs = vs*3 / (rad(j+1)^3-rad(j)^3);
			mu = ro*vs.^2;
			la = ro*vp.^2 - 2.*mu;
			ro = ro*1e3;
			la = la*1e9;
			mu = mu*1e9;

			self.density(j) = ro;
			self.lame_lambda(j) = la;
			self.lame_mu(j) = mu;
		  end

		end % }}}
	end
end

function intnat = naturetointeger(strnat) % {{{
	intnat=zeros(length(strnat),1);
	for i=1:length(strnat),
		switch strnat{i},
		case 'damageice'
			intnat(i)=1;
		case 'estar'
			intnat(i)=2;
		case 'ice'
			intnat(i)=3;
		case 'enhancedice'
			intnat(i)=4;
		%case 'materials' %this case will never happen, kept to ensure equivalent of codes between IoCodeToMaterialsEnum and IoCodeToNatureEnum
		%	intnat(i)=5;
		case 'litho'
			intnat(i)=6;
		case 'hydro'
			intnat(i)=7;
		otherwise
			error('materials constructor error message: nature of the material not supported yet! (''ice'' or ''litho'' or ''hydro'')');
		end
	end
end % }}}
