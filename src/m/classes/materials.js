class materials {
	/**
	 * MATERIALS class definition
	 *
	 * Usage:
	 *     materials = materials();
	 */
	constructor() {//{{{
		this.nature = [];

		let nargs = arguments.length;
		if (nargs == 0) {
			this.nature=['ice'];
		} else {
			this.nature=arguments;
		}

		// Check this is acceptable
		for (let i = 0; i < this.nature.length; ++i) {
			if (!(strcmpi(this.nature[i], 'litho') || strcmpi(this.nature[i], 'ice') || strcmpi(this.nature[i], 'hydro'))) {
				error('materials constructor error message: nature of the material not supported yet! (\'ice\' or \'litho\' or \'hydro\')');
			}
		}

		// Start filling in the dynamic fields (not truly dynamic under JavaScript)
		for (let i = 0; i < length(this.nature); ++i) {
			let nat = this.nature[i];
			switch (nat) {
				case 'ice':
					this.rho_ice = 0;
					this.rho_water = 0;
					this.rho_freshwater = 0;
					this.mu_water = 0;
					this.heatcapacity = 0;
					this.latentheat = 0;
					this.thermalconductivity = 0;
					this.temperateiceconductivity = 0;
					this.effectiveconductivity_averaging = 0;
					this.meltingpoint = 0;
					this.beta = 0;
					this.mixed_layer_capacity = 0;
					this.thermal_exchange_velocity = 0;
					this.rheology_B = 0;
					this.rheology_n = 0;
					this.rheology_law = 0;
					break;
				case 'litho':
					this.numlayers = 0;
					this.radius = 0;
					this.viscosity = 0;
					this.lame_lambda = 0;
					this.lame_mu = 0;
					this.burgers_viscosity = 0;
					this.burgers_mu = 0;
					this.ebm_alpha = 0;
					this.ebm_delta = 0;
					this.ebm_taul = 0;
					this.ebm_tauh = 0;
					this.rheologymodel = 0;
					this.density = 0;
					this.issolid = 0;
					break;
				case 'hydro':
					this.rho_ice = 0;
					this.rho_water = 0;
					this.rho_freshwater = 0;
					break;
				default:
					error('materials constructor error message: nature of the material not supported yet! (\'ice\' or \'litho\' or \'hydro\')');
			}
			this.earth_density = 0;

			// Set default parameters
			this.setdefaultparameters();
		}
	} //}}}
	disp() {//{{{
		console.log(sprintf('   Materials:'));

		for (let i = 0; i < length(this.nature); ++i) {
			let nat = this.nature[i];
			switch (nat) {
				case 'ice':
					console.log(sprintf('   \nIce:'));
					fielddisplay(this,'rho_ice','ice density [kg/m^3]');
					fielddisplay(this,'rho_water','ocean water density [kg/m^3]');
					fielddisplay(this,'rho_freshwater','fresh water density [kg/m^3]');
					fielddisplay(this,'mu_water','water viscosity [N s/m^2]');
					fielddisplay(this,'heatcapacity','heat capacity [J/kg/K]');
					fielddisplay(this,'thermalconductivity','ice thermal conductivity [W/m/K]');
					fielddisplay(this,'temperateiceconductivity','temperate ice thermal conductivity [W/m/K]');
					fielddisplay(this,'meltingpoint','melting point of ice at 1atm in K');
					fielddisplay(this,'latentheat','latent heat of fusion [J/kg]');
					fielddisplay(this,'beta','rate of change of melting point with pressure [K/Pa]');
					fielddisplay(this,'mixed_layer_capacity','mixed layer capacity [W/kg/K]');
					fielddisplay(this,'thermal_exchange_velocity','thermal exchange velocity [m/s]');
					fielddisplay(this,'rheology_B','flow law parameter [Pa s^(1/n)]');
					fielddisplay(this,'rheology_n','Glen\'s flow law exponent');
					fielddisplay(this,'rheology_law','law for the temperature dependance of the rheology: \'None\', \'BuddJacka\', \'Cuffey\', \'CuffeyTemperate\', \'Paterson\', \'Arrhenius\', \'LliboutryDuval\', \'NyeCO2\', or \'NyeH2O\'');
					break;
				case 'litho':
					console.log(sprintf('   \nLitho:'));
					fielddisplay(this,'numlayers','number of layers (default: 2)');
					fielddisplay(this,'radius','array describing the radius for each interface (numlayers+1) [m]');
					fielddisplay(this,'viscosity','array describing each layer\'s viscosity (numlayers) [Pa.s]');
					fielddisplay(this,'lame_lambda','array describing the lame lambda parameter (numlayers) [Pa]');
					fielddisplay(this,'lame_mu','array describing the shear modulus for each layers (numlayers) [Pa]');
					fielddisplay(this,'burgers_viscosity','array describing each layer\'s transient viscosity, only for Burgers rheologies  (numlayers) [Pa.s]');
					fielddisplay(this,'burgers_mu','array describing each layer\'s transient shear modulus, only for Burgers rheologies  (numlayers) [Pa]');

					fielddisplay(this,'ebm_alpha','array describing each layer\'s exponent parameter controlling the shape of shear modulus curve between taul and tauh, only for EBM rheology (numlayers)');
					fielddisplay(this,'ebm_delta','array describing each layer\'s amplitude of the transient relaxation (ratio between elastic rigity to pre-maxwell relaxation rigity), only for EBM rheology (numlayers)');
					fielddisplay(this,'ebm_taul','array describing each layer\'s starting period for transient relaxation, only for EBM rheology  (numlayers) [s]');
					fielddisplay(this,'ebm_tauh','array describing each layer\'s array describing each layer\'s end period for transient relaxation, only for Burgers rheology (numlayers) [s]');

					fielddisplay(this,'rheologymodel','array describing whether we adopt a Maxwell (0), Burgers (1) or EBM (2) rheology (default: 0)');
					fielddisplay(this,'density','array describing each layer\'s density (numlayers) [kg/m^3]');
					fielddisplay(this,'issolid','array describing whether the layer is solid or liquid (default 1) (numlayers)');
					break;
				case 'hydro':
					console.log(sprintf('   \nHydro:'));
					fielddisplay(this,'rho_ice','ice density [kg/m^3]');
					fielddisplay(this,'rho_water','ocean water density [kg/m^3]');
					fielddisplay(this,'earth_density','mantle density [kg/m^3]');
					fielddisplay(this,'rho_freshwater','fresh water density [kg/m^3]');
					break;
				default:
					error('materials constructor error message: nature of the material not supported yet! (\'ice\' or \'litho\' or \'hydro\')');
			}
		}
	} //}}}
	setdefaultparameters() {//{{{
		for (let i = 0; i < this.nature.length; ++i) {
			let nat = this.nature[i];
			switch (nat) {
				case 'ice':
					// Ice density (kg/m^3)
					this.rho_ice = 917;

					// Ocean water density (kg/m^3)
					this.rho_water = 1023

					// Fresh water density (kg/m^3)
					this.rho_freshwater = 1000;

					// Water viscosity (N.s/m^2)
					this.mu_water = 0.001787;

					// Ice heat capacity cp (J/kg/K)
					this.heatcapacity = 2093;

					// Ice latent heat of fusion L (J/kg)
					this.latentheat = 3.34 * 1e5;

					// Ice thermal conductivity (W/m/K)
					this.thermalconductivity = 2.4;

					// Wet ice thermal conductivity (W/m/K)
					this.temperateiceconductivity = 0.24;

					// Computation of effective conductivity
					this.effectiveconductivity_averaging = 1;

					// The melting point of ice at 1 atmosphere of pressure in K
					this.meltingpoint = 273.15;

					// Rate of change of melting point with pressure (K/Pa)
					this.beta = 9.8 * 1e-8;

					// Mixed layer (ice-water interface) heat capacity (J/kg/K)
					this.mixed_layer_capacity = 3974;

					// Thermal exchange velocity (ice-water interface) (m/s)
					this.thermal_exchange_velocity = 1.00 * 1e-4;

					// Rheology law: what is the temperature dependence of B with T
					// available: none, paterson and arrhenius
					this.rheology_law = 'Paterson';

					// Rheology fields default
					this.rheology_B = 1 * 1e8;
					this.rheology_n = 3;
					break;
				case 'litho':
					// We default to a configuration that enables running GIA 
					// solutions using giacaron and/or giaivins
					this.numlayers = 2;

					// Center of the earth (approximation, must not be 0), then the 
					// lab (lithosphere/asthenosphere boundary) then the surface
					// (with 1d3 to avoid numerical singularities)
					this.radius = [1e3, 6278 * 1e3, 6378 * 1e3];

					this.viscosity = [1e21, 1e40]; // Mantle and lithosphere viscosity (respectively) [Pa.s]
					this.lame_mu = [1.45 * 1e11, 6.7 * 1e10]; // (Pa) // Lithosphere and mantle shear modulus (respectively) [Pa]
					this.lame_lambda = this.lame_mu; // (Pa) // Mantle and lithosphere lamba parameter (respectively) [Pa]
					this.burgers_viscosity = [NaN, NaN];
					this.burgers_mu = [NaN, NaN];

					this.ebm_alpha = [NaN, NaN];
					this.ebm_delta = [NaN, NaN];
					this.ebm_taul = [NaN, NaN];
					this.ebm_tauh = [NaN, NaN];
					this.rheologymodel = [0, 0];
					this.density = [5.51 * 1e3, 5.50 * 1e3]; // (Pa) // Mantle and lithosphere density [kg/m^3]
					this.issolid = [1, 1]; // Is layer solid or liquid?
					break;
				case 'hydro':
					// Ice density (kg/m^3)
					this.rho_ice = 917;

					// Ocean water density (kg/m^3)
					this.rho_water = 1023;

					// Fresh water density (kg/m^3)
					this.rho_freshwater = 1000;
					break;
				default:
					error('materials constructor error message: nature of the material not supported yet! (\'ice\' or \'litho\' or \'hydro\')');
			}

			// Average density of the Earth (kg/m^3)
			this.earth_density = 5512;
		}
	} //}}}
	checkconsistency(md, solution, analyses) {//{{{
		for (let i = 0; i < this.nature.length; ++i) {
			let nat = this.nature[i];
			if (nat == 'ice') {
				md = checkfield(md, 'fieldname', 'materials.rho_ice', '>', 0);
				md = checkfield(md, 'fieldname', 'materials.rho_water', '>', 0);
				md = checkfield(md, 'fieldname', 'materials.rho_freshwater', '>', 0);
				md = checkfield(md, 'fieldname', 'materials.mu_water', '>', 0);
				md = checkfield(md, 'fieldname', 'materials.rheology_B', '>', 0, 'timeseries', 1, 'NaN', 1, 'Inf', 1);
				md = checkfield(md, 'fieldname', 'materials.rheology_n', '>', 0, 'size', [md.mesh.numberofelements, 1]);
				md = checkfield(md, 'fieldname', 'materials.rheology_law', 'values', ['None', 'BuddJacka', 'Cuffey', 'CuffeyTemperate', 'Paterson', 'Arrhenius', 'LliboutryDuval', 'NyeCO2', 'NyeH2O']);
			} else if (nat == 'litho') {
				if (!analyses.includes('LoveAnalysis')) {
					return md;
				}

				md = checkfield(md, 'fieldname', 'materials.numlayers', 'NaN', 1, 'Inf', 1, '>', 0, 'numel', 1);
				md = checkfield(md, 'fieldname', 'materials.radius', 'NaN', 1, 'Inf', 1, 'size', [md.materials.numlayers + 1, 1], '>', 0);
				md = checkfield(md, 'fieldname', 'materials.lame_mu', 'NaN', 1, 'Inf', 1, 'size', [md.materials.numlayers, 1], '>=', 0);
				md = checkfield(md, 'fieldname', 'materials.lame_lambda', 'NaN', 1, 'Inf', 1, 'size', [md.materials.numlayers, 1], '>=', 0);
				md = checkfield(md, 'fieldname', 'materials.issolid', 'NaN', 1, 'Inf', 1, 'size', [md.materials.numlayers, 1], '>=', 0, '<', 2);
				md = checkfield(md, 'fieldname', 'materials.density', 'NaN', 1, 'Inf', 1, 'size', [md.materials.numlayers, 1], '>', 0);
				md = checkfield(md, 'fieldname', 'materials.viscosity', 'NaN', 1, 'Inf', 1, 'size', [md.materials.numlayers, 1], '>=', 0);
				md = checkfield(md, 'fieldname', 'materials.rheologymodel', 'NaN', 1, 'Inf', 1, 'size', [md.materials.numlayers, 1], '>=', 0, '<=', 2);
				md = checkfield(md, 'fieldname', 'materials.burgers_viscosity', 'Inf', 1, 'size', [md.materials.numlayers, 1], '>=', 0);
				md = checkfield(md, 'fieldname', 'materials.burgers_mu', 'Inf', 1, 'size', [md.materials.numlayers, 1], '>=', 0);
				md = checkfield(md, 'fieldname', 'materials.ebm_alpha', 'Inf', 1, 'size', [md.materials.numlayers, 1], '>=', 0);
				md = checkfield(md, 'fieldname', 'materials.ebm_delta', 'Inf', 1, 'size', [md.materials.numlayers, 1], '>=', 0);
				md = checkfield(md, 'fieldname', 'materials.ebm_taul', 'Inf', 1, 'size', [md.materials.numlayers, 1], '>=', 0);
				md = checkfield(md, 'fieldname', 'materials.ebm_tauh', 'Inf', 1, 'size', [md.materials.numlayers, 1], '>=', 0);

				for (let i = 0; i < md.materials.numlayers; ++i) {
					if (md.materials.rheologymodel[i] == 1 && (isnan(md.materials.burgers_viscosity[i]) || isnan(md.materials.burgers_mu[i]))) {
						error('materials checkconsistency error message: Litho burgers_viscosity or burgers_mu has NaN values, inconsistent with rheologymodel choice');
					}
					if (md.materials.rheologymodel[i] == 2 && (isnan(md.materials.ebm_alpha[i]) || isnan(md.materials.ebm_delta[i])) || isnan(md.materials.ebm_taul[i] || isnan(md.materials.ebm_tauh[i]))) {
						error('materials checkconsistency error message: Litho ebm_alpha, ebm_delta, ebm_taul or ebm_tauh has NaN values, inconsistent with rheologymodel choice');
					}
				}
				if (md.materials.issolid[0] == 0 || md.materials.lame_mu[0] == 0) {
					error('First layer must be solid (issolid[0] > 0 AND lame_mu[0] > 0). Add a weak inner core if necessary.');
				}
				let ind = find(md.materials.issolid == 0);
				if (sum(ismember(diff(ind), 1) >= 1)) { // If there are at least two consecutive indices that contain issolid = 0
					error('Fluid layers detected at layers #' + ind + ', but having 2 or more adjacent fluid layers is not supported yet. Consider merging them.');
				}
			} else if (nat == 'hydro') {
				md = checkfield(md, 'fieldname', 'materials.rho_ice', '>', 0);
				md = checkfield(md, 'fieldname', 'materials.rho_water', '>', 0);
				md = checkfield(md, 'fieldname', 'materials.earth_density', '>', 0, 'numel', 1);
				md = checkfield(md, 'fieldname', 'materials.rho_freshwater', '>', 0);
			} else {
				error('materials constructor error message: nature of the material not supported yet! (\'ice\' or \'litho\' or \'hydro\')');
			}
		}

		return md;
	} //}}}
	marshall(md, prefix, fid) {//{{{
		// 1: MatdamageiceEnum 2: MatestarEnum 3: MaticeEnum 4: MatenhancediceEnum 5: MaterialsEnum
		WriteData(fid, prefix, 'name', 'md.materials.nature', 'data', naturetointeger(this.nature), 'format', 'IntMat', 'mattype', 3);
		WriteData(fid, prefix, 'name', 'md.materials.type', 'data', 5, 'format', 'Integer'); // DANGER: this can evolve if you have classes
		for (let i = 0; i < this.nature.length; ++i) {
			let nat = this.nature[i];
			if (nat == 'ice') {
				WriteData(fid, prefix, 'object', this, 'class', 'materials', 'fieldname', 'rho_ice', 'format', 'Double');
				WriteData(fid, prefix, 'object', this, 'class', 'materials', 'fieldname', 'rho_water', 'format', 'Double');
				WriteData(fid, prefix, 'object', this, 'class', 'materials', 'fieldname', 'rho_freshwater', 'format', 'Double');
				WriteData(fid, prefix, 'object', this, 'class', 'materials', 'fieldname', 'mu_water', 'format', 'Double');
				WriteData(fid, prefix, 'object', this, 'class', 'materials', 'fieldname', 'heatcapacity', 'format', 'Double');
				WriteData(fid, prefix, 'object', this, 'class', 'materials', 'fieldname', 'latentheat', 'format', 'Double');
				WriteData(fid, prefix, 'object', this, 'class', 'materials', 'fieldname', 'thermalconductivity', 'format', 'Double');
				WriteData(fid, prefix, 'object', this, 'class', 'materials', 'fieldname', 'temperateiceconductivity', 'format', 'Double');
				WriteData(fid, prefix, 'object', this, 'class', 'materials', 'fieldname', 'effectiveconductivity_averaging', 'format', 'Integer');
				WriteData(fid, prefix, 'object', this, 'class', 'materials', 'fieldname', 'meltingpoint', 'format', 'Double');
				WriteData(fid, prefix, 'object', this, 'class', 'materials', 'fieldname', 'beta', 'format', 'Double');
				WriteData(fid, prefix, 'object', this, 'class', 'materials', 'fieldname', 'mixed_layer_capacity', 'format', 'Double');
				WriteData(fid, prefix, 'object', this, 'class', 'materials', 'fieldname', 'thermal_exchange_velocity', 'format', 'Double');
				WriteData(fid, prefix, 'object', this, 'class', 'materials', 'fieldname', 'rheology_B', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts);
				WriteData(fid, prefix, 'object', this, 'class', 'materials', 'fieldname', 'rheology_n', 'format', 'DoubleMat', 'mattype', 2);
				WriteData(fid, prefix, 'data', this.rheology_law, 'name', 'md.materials.rheology_law', 'format', 'String');
			} else if (nat == 'litho') {
				WriteData(fid, prefix, 'object', this, 'class', 'materials', 'fieldname', 'numlayers', 'format', 'Integer');
				WriteData(fid, prefix, 'object', this, 'class', 'materials', 'fieldname', 'radius', 'format', 'DoubleMat', 'mattype', 3);
				WriteData(fid, prefix, 'object', this, 'class', 'materials', 'fieldname', 'lame_mu', 'format', 'DoubleMat', 'mattype', 3);
				WriteData(fid, prefix, 'object', this, 'class', 'materials', 'fieldname', 'lame_lambda', 'format', 'DoubleMat', 'mattype', 3);
				WriteData(fid, prefix, 'object', this, 'class', 'materials', 'fieldname', 'issolid', 'format', 'DoubleMat', 'mattype', 3);
				WriteData(fid, prefix, 'object', this, 'class', 'materials', 'fieldname', 'density', 'format', 'DoubleMat', 'mattype', 3);
				WriteData(fid, prefix, 'object', this, 'class', 'materials', 'fieldname', 'viscosity', 'format', 'DoubleMat', 'mattype', 3);
				WriteData(fid, prefix, 'object', this, 'class', 'materials', 'fieldname', 'rheologymodel', 'format', 'DoubleMat', 'mattype', 3);
				WriteData(fid, prefix, 'object', this, 'class', 'materials', 'fieldname', 'burgers_viscosity', 'format', 'DoubleMat', 'mattype', 3);
				WriteData(fid, prefix, 'object', this, 'class', 'materials', 'fieldname', 'burgers_mu', 'format', 'DoubleMat', 'mattype', 3);
				WriteData(fid, prefix, 'object', this, 'class', 'materials', 'fieldname', 'ebm_alpha', 'format', 'DoubleMat', 'mattype', 3);
				WriteData(fid, prefix, 'object', this, 'class', 'materials', 'fieldname', 'ebm_delta', 'format', 'DoubleMat', 'mattype', 3);
				WriteData(fid, prefix, 'object', this, 'class', 'materials', 'fieldname', 'ebm_taul', 'format', 'DoubleMat', 'mattype', 3);
				WriteData(fid, prefix, 'object', this, 'class', 'materials', 'fieldname', 'ebm_tauh', 'format', 'DoubleMat', 'mattype', 3);
				// Compute earth density compatible with our layer density distribution
				let earth_density = 0;
				for (let i = 0; i < this.numlayers; ++i) {
					earth_density = earth_density + (Math.pow(this.radius[i + 1], 3) - Math.pow(this.radius[i], 3)) * this.density[i];
				}
				earth_density = earth_density / Math.pow(this.radius[this.numlayers]);
				this.earth_density = earth_density;
			} else if (nat == 'hydro') {
				WriteData(fid, prefix, 'object', this, 'class', 'materials', 'fieldname', 'rho_ice', 'format', 'Double');
				WriteData(fid, prefix, 'object', this, 'class', 'materials', 'fieldname', 'rho_water', 'format', 'Double');
				WriteData(fid, prefix, 'object', this, 'class', 'materials', 'fieldname', 'rho_freshwater', 'format', 'Double');
			} else {
				error('materials constructor error message: nature of the material not supported yet! (\'ice\' or \'litho\' or \'hydro\')');
			}
			WriteData(fid, prefix, 'data', this.earth_density, 'name', 'md.materials.earth_density', 'format', 'Double');
		}
	} //}}}
	extrude(md) {//{{{
		for (let i = 0; i < this.nature.length; ++i) {
			let nat = this.nature[i];
			if (nat == 'ice') {
				this.rheology_B = project3d(md, 'vector', this.rheology_B, 'type', 'node');
				this.rheology_n = project3d(md, 'vector', this.rheology_n, 'type', 'element');
			}
		}
	} //}}}
} //}}}

function naturetointeger(strnat) {//{{{
	let intnat = zeros(strnat.length, 1);
	for (let i = 0; i < strnat.length; ++i) {
		switch (strnat[i]) {
			case 'damageice':
				intnat[i] = 1;
				break;
			case 'estar':
				intnat[i] = 2;
				break;
			case 'ice':
				intnat[i] = 3;
				break;
			case 'enhancedice':
				intnat[i] = 4;
				break;
			//case 'materials': // This case will never happen, kept to ensure equivalent of codes between IoCodeToMaterialsEnum and IoCodeToNatureEnum
			//	intnat[i] = 5;
			//	break;
			case 'litho':
				intnat[i] = 6;
				break;
			case 'hydro':
				intnat[i] = 7;
				break;
			default:
				error('materials constructor error message: nature of the material not supported yet! (\'ice\' or \'litho\' or \'hydro\')');
		}
	}
	return intnat;
}// }}}
