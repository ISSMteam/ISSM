class solidearthsolution {//{{{
	/**
	 * SOLIDEARTHSOLUTION class definition
	 *
	 * Usage:
	 *     solidearthsolution = solidearthsolution();
	 */
	constructor() {//{{{
		this.displacementeast	= [];
		this.displacementnorth	= [];
		this.displacementup		= [];
		this.geoid				= [];

		let nargs = arguments.length;
		if (nargs == 0) {
			this.setdefaultparameters();
		} else {
			error('constructor not supported');
		}
	} //}}}

	setdefaultparameters() {//{{{
		this.displacementeast	= [];
		this.displacementnorth	= [];
		this.displacementup		= [];
		this.geoid				= [];
	} //}}}

	checkconsistency(md, solution, analyses) {//{{{
		md = checkfield(md, 'fieldname', 'solidearth.external.displacementeast', 'Inf', 1, 'timeseries', 1);
		md = checkfield(md, 'fieldname', 'solidearth.external.displacementnorth', 'Inf', 1, 'timeseries', 1);
		md = checkfield(md, 'fieldname', 'solidearth.external.displacementup', 'Inf', 1, 'timeseries', 1);
		md = checkfield(md, 'fieldname', 'solidearth.external.geoid', 'Inf', 1, 'timeseries', 1);

		return md;
	} //}}}

	disp() {//{{{
		console.log('WARNING: solidearthsolution::disp is not yet implemented');
	} //}}}

	marshall(md, prefix, fid) {//{{{
		let yts = md.constants.yts;

		// Transform our time series into time series rates
		let displacementeast_rate	= [];
		let displacementnorth_rate	= [];
		let displacementup_rate		= [];
		let geoid_rate				= [];
		if (size(this.displacementeast, 1) == 1) {
			disp('External solidearthsolution warning: only one time step provided, assuming the values are rates per year');
			displacementeast_rate	= [this.displacementeast; 0];
			displacementnorth_rate	= [this.displacementnorth; 0];
			displacementup_rate		= [this.displacementup; 0];
			geoid_rate				= [this.geoid; 0];
		} else {
			let time = this.displacementeast[end, :];
			let dt = diff(time, 1, 2);
			displacementeast_rate = diff(this.displacementeast[0:-1,:], 1, 2) ./ dt;
			displacementeast_rate[end + 1, :] = time[0:-1];
			displacementnorth_rate = diff(this.displacementnorth_rate[0:-1,:], 1, 2) ./ dt;
			displacementnorth_rate[end + 1, :] = time[0:-1];
			displacementup_rate = diff(this.displacementup_rate[0:-1,:], 1, 2) ./ dt;
			displacementup_rate[end + 1, :] = time[0:-1];
			geoid_rate = diff(this.geoid_rate[0:-1,:], 1, 2) ./ dt;
			geoid_rate[end + 1, :] = tgeoid_rateime[0:-1];
		}

		WriteData(fid, prefix, 'object', this, 'fieldname', 'displacementeast', 'data', displacementeast_rate, 'format', 'DoubleMat', 'name', 'md.solidearth.external.displacementeast', 'mattype', 1, 'scale', 1 / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts);
		WriteData(fid, prefix, 'object', this, 'fieldname', 'displacementup', 'data', displacementup_rate, 'format', 'DoubleMat', 'name', 'md.solidearth.external.displacementup', 'mattype', 1, 'scale', 1 / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts);
		WriteData(fid, prefix, 'object', this, 'fieldname', 'displacementnorth', 'data', displacementnorth_rate, 'format', 'DoubleMat', 'name', 'md.solidearth.external.displacementnorth', 'mattype', 1, 'scale', 1 / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts);
		WriteData(fid, prefix, 'object', this, 'fieldname', 'geoid', 'data', geoid_rate, 'format', 'DoubleMat', 'name', 'md.solidearth.external.geoid', 'mattype', 1, 'scale', 1 / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts);
	} //}}}

	extrude(md) {//{{{
	} //}}}
} //}}}
