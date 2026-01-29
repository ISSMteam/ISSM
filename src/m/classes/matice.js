//MATICE class definition
//
//   Usage:
//      matice = new matice();

function matice(){
	//methods
	this.setdefaultparameters = function(){ // {{{

		//ice density (kg/m^3)
		this.rho_ice=917;

		//ocean water density (kg/m^3)
		this.rho_water=1023;

		//fresh water density (kg/m^3)
		this.rho_freshwater=1000;

		//water viscosity (N.s/m^2)
		this.mu_water=0.001787;

		//ice heat capacity cp (J/kg/K)
		this.heatcapacity=2093;

		//ice latent heat of fusion L (J/kg)
		this.latentheat=3.34e5;

		//ice thermal conductivity (W/m/K)
		this.thermalconductivity=2.4;

		//wet ice thermal conductivity (W/m/K)
		this.temperateiceconductivity=.24;

		//computation of effective conductivity
		this.effectiveconductivity_averaging=1;

		//the melting point of ice at 1 atmosphere of pressure in K
		this.meltingpoint=273.15;

		//rate of change of melting point with pressure (K/Pa)
		this.beta=9.8e-8;

		//mixed layer (ice-water interface) heat capacity (J/kg/K)
		this.mixed_layer_capacity=3974;

		//thermal exchange velocity (ice-water interface) (m/s)
		this.thermal_exchange_velocity=1e-4;

		//Rheology law: what is the temperature dependence of B with T
		//available: none, paterson and arrhenius
		this.rheology_law='Paterson';

		// Rheology for ice
		this.rheology_B = 2.1e8
		this.rheology_n = 3

		//SLC
		this.earth_density= 5512; // average density of the Earth, (kg/m^3)


	} //}}}
	this.disp = function() {// {{{
		console.log(sprintf('   Materials:'));

		fielddisplay(this,'rho_ice','ice density [kg/m^3]');
		fielddisplay(this,'rho_water','ocean water density [kg/m^3]');
		fielddisplay(this,'rho_freshwater','fresh water density [kg/m^3]');
		fielddisplay(this,'mu_water','water viscosity [N s/m^2]');
		fielddisplay(this,'heatcapacity','heat capacity [J/kg/K]');
		fielddisplay(this,'thermalconductivity','ice thermal conductivity [W/m/K]');
		fielddisplay(this,'temperateiceconductivity','temperate ice thermal conductivity [W/m/K]');
		fielddisplay(this,'effectiveconductivity_averaging','computation of effectiveconductivity: (0) arithmetic mean, (1) harmonic mean, (2) geometric mean (default)');
		fielddisplay(this,'meltingpoint','melting point of ice at 1atm in K');
		fielddisplay(this,'latentheat','latent heat of fusion [J/kg]');
		fielddisplay(this,'beta','rate of change of melting point with pressure [K/Pa]');
		fielddisplay(this,'mixed_layer_capacity','mixed layer capacity [W/kg/K]');
		fielddisplay(this,'thermal_exchange_velocity','thermal exchange velocity [m/s]');
		fielddisplay(this,'rheology_B','flow law parameter [Pa s^(1/n)]');
		fielddisplay(this,'rheology_n','Glen\'s flow law exponent');
		fielddisplay(this,'rheology_law','law for the temperature dependance of the rheology: \'None\', \'BuddJacka\', \'Cuffey\', \'CuffeyTemperate\', \'Paterson\', \'Arrhenius\', \'LliboutryDuval\', \'NyeH2O\', or \'NyeCO2\'');
		fielddisplay(this,'earth_density','Mantle density [kg/m^-3]');

	} // }}}
	this.extrude = function(md) {//{{{
		this.rheology_B=project3d(md,'vector',this.rheology_B,'type','node');
		this.rheology_n=project3d(md,'vector',this.rheology_n,'type','element');
		return this;
	}//}}}
	this.classname = function() {// {{{
		return 'matice';
	} // }}}
	this.checkconsistency = function(md,solution,analyses) { // {{{
		if (solution == 'TransientSolution' && md.transient.isslc) {
			checkfield(md,'fieldname','materials.earth_density','>',0,'numel',[1])
		} else {
			checkfield(md,'fieldname','materials.rho_ice','>',0);
			checkfield(md,'fieldname','materials.rho_water','>',0);
			checkfield(md,'fieldname','materials.rho_freshwater','>',0);
			checkfield(md,'fieldname','materials.mu_water','>',0);
			checkfield(md,'fieldname','materials.rheology_B','>',0,'size','universal','NaN',1,'Inf',1);
			checkfield(md,'fieldname','materials.rheology_n','>',0,'size','universal','NaN',1,'Inf',1);
			checkfield(md,'fieldname','materials.rheology_law','values',['None','BuddJacka','Cuffey','CuffeyTemperate','Paterson','Arrhenius','LliboutryDuval','NyeH2O','NyeCO2']);
			checkfield(md,'fieldname','materials.effectiveconductivity_averaging','numel',[1],'values',[0,1,2])
		}
	} // }}}
	this.marshall=function(md,prefix,fid) { //{{{
		WriteData(fid,prefix,'name','md.materials.type','data',3,'format','Integer');
		WriteData(fid,prefix,'object',this,'class','materials','fieldname','rho_ice','format','Double');
		WriteData(fid,prefix,'object',this,'class','materials','fieldname','rho_water','format','Double');
		WriteData(fid,prefix,'object',this,'class','materials','fieldname','rho_freshwater','format','Double');
		WriteData(fid,prefix,'object',this,'class','materials','fieldname','mu_water','format','Double');
		WriteData(fid,prefix,'object',this,'class','materials','fieldname','heatcapacity','format','Double');
		WriteData(fid,prefix,'object',this,'class','materials','fieldname','latentheat','format','Double');
		WriteData(fid,prefix,'object',this,'class','materials','fieldname','thermalconductivity','format','Double');
		WriteData(fid,prefix,'object',this,'class','materials','fieldname','temperateiceconductivity','format','Double');
		WriteData(fid,prefix,'object',this,'class','materials','fieldname','effectiveconductivity_averaging','format','Integer');
		WriteData(fid,prefix,'object',this,'class','materials','fieldname','meltingpoint','format','Double');
		WriteData(fid,prefix,'object',this,'class','materials','fieldname','beta','format','Double');
		WriteData(fid,prefix,'object',this,'class','materials','fieldname','mixed_layer_capacity','format','Double');
		WriteData(fid,prefix,'object',this,'class','materials','fieldname','thermal_exchange_velocity','format','Double');
		if ((shape(self.rheology_B)[0] == md.mesh.numberofvertices) || (shape(self.rheology_B)[0] == md.mesh.numberofvertices + 1) || (shape(self.rheology_B)[0] == md.mesh.numberofelements && shape(self.rheology_B)[1] > 1)) {
			mattype = 1
			tsl = md.mesh.numberofvertices
		} else {
			mattype = 2
			tsl = md.mesh.numberofelements
		}
		WriteData(fid,prefix,'object',this,'class','materials','fieldname','rheology_B','format','DoubleMat','mattype',mattype,'timeserieslength',tsl+1,'yts',md.constants.yts);
		WriteData(fid,prefix,'object',this,'class','materials','fieldname','rheology_n','format','DoubleMat','mattype',2);
		WriteData(fid,prefix,'data',this.rheology_law,'name','md.materials.rheology_law','format','String');
		WriteData(fid,prefix,'object',this,'class','materials','fieldname','earth_density','format','Double');

	}//}}}
	this.fix=function() { //{{{
	}//}}}
	//properties
	// {{{
	this.rho_ice							= 0;
	this.rho_water							= 0;
	this.rho_freshwater						= 0;
	this.mu_water							= 0;
	this.heatcapacity						= 0;
	this.latentheat							= 0;
	this.thermalconductivity				= 0;
	this.temperateiceconductivity			= 0;
	this.effectiveconductivity_averaging	= 0;
	this.meltingpoint						= 0;
	this.beta								= 0;
	this.mixed_layer_capacity				= 0;
	this.thermal_exchange_velocity			= 0;
	this.rheology_B							= NaN;
	this.rheology_n							= NaN;
	this.rheology_law						= '';

	//SLC
	this.earth_density						= 0;

	this.setdefaultparameters();
	//}}}
}
