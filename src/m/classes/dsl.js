class dsl {//{{{
	/**
	 * DSL class definition
	 *
	 * Usage:
	 *     dsl=dsl();
	 * 
	 * TODO:
	 * - Implement from struct constructor (see dsl.m)?
	 */
	constructor() {//{{{
		this.global_average_thermosteric_sea_level; // Corresponds to zostoga field in CMIP5 archives. Specified as a temporally variable quantity (in m)
		this.sea_surface_height_above_geoid; // Corresponds to zos field in CMIP5 archives. Spatial average is 0. Specified as a spatio-temporally variable quantity (in m)
		this.sea_water_pressure_at_sea_floor; // Corresponds to bpo field in CMIP5 archives. Specified as a spatio-temporally variable quantity (in m equivalent, not in Pa!)

		let nargs = arguments.length;
		if (nargs == 0) {
			this.setdefaultparameters();
		} else {
			error('constructor not supported');
		}
	} //}}}

	setdefaultparameters() {//{{{
		this.global_average_thermosteric_sea_level = NaN;
		this.sea_surface_height_above_geoid = NaN;
		this.sea_water_pressure_at_sea_floor = NaN;
	} //}}}

	checkconsistency(md, solution, analyses) {//{{{
		if (!analyses.includes('SealevelchangeAnalysis') || (solution === 'TransientSolution' && !md.transient.isslc) || !md.transient.isoceantransport) {
			return md;
		}

		md = checkfield(md, 'fieldname', 'dsl.global_average_thermosteric_sea_level', 'NaN', 1, 'Inf', 1);
		md = checkfield(md, 'fieldname', 'dsl.sea_surface_height_above_geoid', 'NaN', 1, 'Inf', 1, 'timeseries', 1);
		md = checkfield(md, 'fieldname', 'dsl.sea_water_pressure_at_sea_floor', 'NaN', 1, 'Inf', 1, 'timeseries', 1);

		if (md.solidearth.settings.compute_bp_grd) {
			md = checkfield(md, 'fieldname', 'dsl.sea_water_pressure_at_sea_floor', 'empty', 1);
		}

		return md;
	} //}}}

	disp() {//{{{
		console.log('WARNING: dsl::disp is not yet implemented');
	} //}}}

	marshall(md, prefix, fid) {//{{{
		WriteData(fid, prefix, 'name', 'md.dsl.model', 'data', 1, 'format', 'Integer');
		WriteData(fid, prefix, 'object', this, 'fieldname', 'global_average_thermosteric_sea_level', 'format', 'DoubleMat', 'mattype', 2, 'timeseries', 1, 'timeserieslength', 2, 'yts', md.constants.yts); // mattype 2, because we are sending a GMSL value identical everywhere on each element
		WriteData(fid, prefix, 'object', this, 'fieldname', 'sea_surface_height_above_geoid', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts); // mattype 1 because we specify DSL at vertex locations
		WriteData(fid, prefix, 'object', this, 'fieldname', 'sea_water_pressure_at_sea_floor', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts); // mattype 1 because we specify bottom pressure at vertex locations
	} //}}}

	extrude(md) {//{{{
		this.sea_surface_height_above_geoid = project3d(md, 'vector', this.sea_surface_height_above_geoid, 'type', 'node', 'layer', 1);
		this.sea_water_pressure_at_sea_floor = project3d(md, 'vector', this.sea_water_pressure_at_sea_floor, 'type', 'node', 'layer', 1);
	} //}}}
} //}}}
