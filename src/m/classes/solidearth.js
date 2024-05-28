class solidearth {//{{{
	/**
	 * SOLIDEARTH class definition
	 *
	 * Usage:
	 *     solidearth = solidearth();
	 *     solidearth = solidearth('earth');
	 */
	constructor() {//{{{
		this.settings = new solidearthsettings();
		this.external = null;
		this.lovenumbers = new lovenumbers();
		this.rotational = new rotational();
		this.planetradius = planetradius('earth');
		this.requested_outputs = [];
		this.transitions = [];
		this.partitionice = [];
		this.partitionhydro = [];
		this.partitionocean = [];

		let nargs = arguments.length;
		if (nargs == 0) {
			this.setdefaultparameters('earth');
		} else if (nargs == 1) {
			this.setdefaultparameters(arguments[0]);
		} else {
			error('solidearth constructor error message: zero or one argument only!');
		}
	} //}}}

	disp() {//{{{
		console.log('WARNING: solidearth::disp is not yet implemented');
	} //}}}

	setdefaultparameters(planet) {//{{{
		// Output default
		this.requested_outputs = ['default'];

		// Transitions should be an array
		this.transitions = [];

		// No partitions requested for barystatic contribution
		this.partitionice = [];
		this.partitionhydro = [];
		this.partitionocean = [];

		// No external solutions by default
		this.external = null;

		// Planet radius
		this.planetradius = planetradius(planet);
	} //}}}

	checkconsistency(md, solution, analyses) {//{{{
		if (!analyses.includes('SealevelchangeAnalysis') || (solution === 'TransientSolution' && !md.transient.isslc)) {
			return md;
		}

		md = checkfield(md, 'fieldname', 'solidearth.requested_outputs', 'stringrow', 1);

		this.settings.checkconsistency(md, solution, analyses);
		this.lovenumbers.checkconsistency(md, solution, analyses);
		this.rotational.checkconsistency(md, solution, analyses);

		if (this.external != null) {
			if (typeof(this.external) != 'solidearthsolution') {
				error('solidearth consistency check: external field should be a solidearthsolution');
			}
		}

		return md;
	} //}}}

	defaultoutputs(md) {//{{{
		return ['Sealevel'];
	} //}}}

	marshall(md, prefix, fid) {//{{{
		WriteData(fid, prefix, 'object', this, 'fieldname', 'planetradius', 'format', 'Double');
		WriteData(fid, prefix, 'object', this, 'fieldname', 'transitions', 'format', 'MatArray');

		let npartice = 0;
		if (this.partitionice.length) {
			npartice = max(this.partitionice) + 2;
		}

		let nparthydro = 0;
		if (this.partitionhydro.length) {
			nparthydro = max(this.partitionhydro) + 2;
		}

		let npartocean = 0;
		if (this.partitionocean.length) {
			npartocean = max(this.partitionocean) + 2;
		}

		WriteData(fid, prefix, 'object', this, 'fieldname', 'partitionice', 'mattype', 1, 'format', 'DoubleMat');
		WriteData(fid, prefix, 'data', npartice, 'format', 'Integer', 'name', 'md.solidearth.npartice');
		WriteData(fid, prefix, 'object', this, 'fieldname', 'partitionhydro', 'mattype', 1, 'format', 'DoubleMat');
		WriteData(fid, prefix, 'data', nparthydro, 'format', 'Integer', 'name', 'md.solidearth.nparthydro');
		WriteData(fid, prefix, 'object', this, 'fieldname', 'partitionocean', 'mattype', 1, 'format', 'DoubleMat');
		WriteData(fid, prefix, 'data', npartocean, 'format', 'Integer', 'name', 'md.solidearth.npartocean');

		this.settings.marshall(md, prefix, fid);
		this.lovenumbers.marshall(md, prefix, fid);
		this.rotational.marshall(md, prefix, fid);

		if (this.external != null) {
			WriteData(fid, prefix, 'data', 1, 'format', 'Integer', 'name', 'md.solidearth.isexternal');
			this.external.marshall(md, prefix, fid);
		} else {
			WriteData(fid, prefix, 'data', 0, 'format', 'Integer', 'name', 'md.solidearth.isexternal');
		}

		// Process requested outputs
		let outputs = this.requested_outputs;
		let pos = find(ismember(outputs, 'default'));
		if (pos.length) {
			/*
			NOTE: In order to handle case where user has added 'default' more 
			than once to this.requested_outputs, need to remove elements by 
			index in reverse order.
			*/
			for (let i = (pos.length - 1); i >= 0; --i) {
				outputs.splice(pos[i], 1);
			}

			outputs.push(this.defaultoutputs(md));
		}
		WriteData(fid, prefix, 'data', outputs, 'name', 'md.solidearth.requested_outputs', 'format', 'StringArray');
	} //}}}

	extrude(md) {//{{{
	} //}}}
} //}}}
