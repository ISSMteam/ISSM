class hydrologyshreve {//{{{
	/**
	 * HYDROLOGYSHREVE class definition
	 *
	 * Usage:
	 *     hydrologyshreve = hydrologyshreve();
	 */
	constructor() {//{{{
		this.spcwatercolumn = NaN;
		this.stabilization = 0;
		this.requested_outputs = [];

		if (arguments.length == 0) {
			this.setdefaultparameters();
		} else {
			error('constructor not supported');
		}
	} //}}}

	disp() {//{{{
		console.log(sprintf('   hydrologyshreve solution parameters:'));
		fielddisplay(this, 'spcwatercolumn', 'water thickness constraints (NaN means no constraint) [m]');
		fielddisplay(this, 'stabilization', 'artificial diffusivity (default: 1). can be more than 1 to increase diffusivity.');
		fielddisplay(this, 'requested_outputs', 'additional outputs requested');
	} //}}}

	setdefaultparameters() {//{{{
		// Type of stabilization to use 0:nothing 1:artificial_diffusivity
		this.stabilization = 1;
		this.requested_outputs = ['default'];
	} //}}}

	checkconsistency(md, solution, analyses) {//{{{
		// Early return
		if (!analyses.contains('HydrologyShreveAnalysis') || (solution == 'TransientSolution' && !md.transient.ishydrology)) {
			return md;
		}
		md = checkfield(md, 'fieldname', 'hydrology.spcwatercolumn', 'Inf', 1, 'timeseries', 1);
		md = checkfield(md, 'fieldname', 'hydrology.stabilization', '>=', 0);
		return md;
	} //}}}

	defaultoutputs(md) {//{{{
		return ['Watercolumn', 'HydrologyWaterVx', 'HydrologyWaterVy'];
	} //}}}

	marshall(md, prefix, fid) {//{{{
		WriteData(fid, prefix, 'name', 'md.hydrology.model', 'data', 2, 'format', 'Integer');
		WriteData(fid, prefix, 'object', this, 'fieldname', 'spcwatercolumn', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts);
		WriteData(fid, prefix, 'object', this, 'fieldname', 'stabilization', 'format', 'Double');

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
			let default_outputs = this.defaultoutputs(md);
			for (let i = 0; i < default_outputs.length; ++i) {
				outputs.push(default_outputs[i]);
			}
		}
		WriteData(fid, prefix, 'data', outputs, 'name', 'md.hydrology.requested_outputs', 'format', 'StringArray');
	} //}}}

	extrude(md) {//{{{
	} //}}}
} //}}}
