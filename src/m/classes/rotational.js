class rotational {//{{{
	/**
	 * ROTATIONAL class definition
	 *
	 * Usage:
	 *     rotational = rotational();
	 */
	constructor() {//{{{
		this.equatorialmoi		= 0;
		this.polarmoi			= 0;
		this.angularvelocity	= 0;

		let nargs = arguments.length;
		if (nargs == 0) {
			this.setdefaultparameters();
		} else {
			error('constructor not supported');
		}
	} //}}}

	disp() {//{{{
		console.log('WARNING: rotational::disp is not yet implemented');
	} //}}}

	setdefaultparameters() {//{{{
		// Moment of inertia
		this.equatorialmoi	= 8.0077e37; // [kg m^2]
		this.polarmoi		= 8.0345e37; // [kg m^2]

		// Mean rotational velocity of earth
		this.angularvelocity = 7.2921e-5; // [s^-1]
	} //}}}

	checkconsistency(md, solution, analyses) {//{{{
		if (!analyses.includes('SealevelchangeAnalysis') || (solution === 'TransientSolution' && !md.transient.isslc)) {
			return md;
		}

		md = checkfield(md, 'fieldname', 'solidearth.rotational.equatorialmoi', 'NaN', 1, 'Inf', 1);
		md = checkfield(md, 'fieldname', 'solidearth.rotational.polarmoi', 'NaN', 1, 'Inf', 1);
		md = checkfield(md, 'fieldname', 'solidearth.rotational.angularvelocity', 'NaN', 1, 'Inf', 1);

		return md;
	} //}}}

	defaultoutputs(md) {//{{{
		return [];
	} //}}}

	marshall(md, prefix, fid) {//{{{
		WriteData(fid, prefix, 'object', this, 'fieldname', 'equatorialmoi', 'name', 'md.solidearth.rotational.equatorialmoi', 'format', 'Double');
		WriteData(fid, prefix, 'object', this, 'fieldname', 'polarmoi', 'name', 'md.solidearth.rotational.polarmoi', 'format', 'Double');
		WriteData(fid, prefix, 'object', this, 'fieldname', 'angularvelocity', 'name', 'md.solidearth.rotational.angularvelocity', 'format', 'Double');
	} //}}}

	extrude(md) {//{{{
	} //}}}
} //}}}
