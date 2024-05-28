class masstransport {//{{{
	/**
	 * MASSTRANSPORT class definition
	 *
	 * Usage:
	 *     masstransport = masstransport();
	 */
	constructor() {//{{{
		this.spcthickness = NaN;
		this.isfreesurface = 0;
		this.min_thickness = 0;
		this.hydrostatic_adjustment = 0;
		this.stabilization = 0;
		this.vertex_pairing = 0;
		this.penalty_factor = 0;
		this.requested_outputs = 0;

		if (arguments.length == 0) {
			this.setdefaultparameters();
		} else {
			error('constructor not supported');
		}
	} //}}}

	disp() {//{{{
		console.log('WARNING: masstransport::disp is not yet implemented');
	} //}}}

	setdefaultparameters() {//{{{
		// Type of stabilization to use 0:nothing 1:artificial_diffusivity 3:Discontinuous Galerkin
		this.stabilization = 1;

		// Factor applied to compute the penalties kappa=max(stiffness matrix)*10^penalty_factor
		this.penalty_factor = 3;

		// Minimum ice thickness that can be used
		this.min_thickness = 1;

		// Hydrostatic adjustment
		this.hydrostatic_adjustment = 'Absolute';

		// Default output
		this.requested_outputs = ['default'];
	} //}}}

	checkconsistency(md, solution, analyses) {//{{{
		// Early return
		if (analyses.includes('HydrologyShreveAnalysis') || (solution == 'TransientSolution' && !md.trans.ismasstransport)) {
			return md;
		}

		md = checkfield(md, 'fieldname', 'masstransport.spcthickness', 'Inf', 1, 'timeseries', 1);
		md = checkfield(md, 'fieldname', 'masstransport.isfreesurface', 'values', [0, 1]);
		md = checkfield(md, 'fieldname', 'masstransport.hydrostatic_adjustment', 'values', ['Absolute', 'Incremental']);
		md = checkfield(md, 'fieldname', 'masstransport.stabilization', 'values', [0, 1, 2, 3, 4, 5]);
		md = checkfield(md, 'fieldname', 'masstransport.min_thickness', '>', 0);
		md = checkfield(md, 'fieldname', 'masstransport.requested_outputs', 'stringrow', 1);
		if (!any(isnan(md.stressbalance.vertex.vertex_pairing))) {
			md = checkfield(md, 'fieldname', 'stressbalance.vertex_pairing', '>', 0);
		}

		return md;
	} //}}}

	defaultoutputs(md) {//{{{
		return ['Thickness', 'Surface', 'Base'];
	} //}}}

	marshall(md, prefix, fid) {//{{{
		let yts = md.constants.yts;

		WriteData(fid, prefix, 'object', this, 'fieldname', 'spcthickness', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts);
		WriteData(fid, prefix, 'object', this, 'fieldname', 'isfreesurface', 'format', 'Boolean');
		WriteData(fid, prefix, 'object', this, 'fieldname', 'min_thickness', 'format', 'Double');
		WriteData(fid, prefix, 'data', this.hydrostatic_adjustment, 'format', 'String', 'name', 'md.masstransport.hydrostatic_adjustment');
		WriteData(fid, prefix, 'object', this, 'fieldname', 'stabilization', 'format', 'Integer');
		WriteData(fid, prefix, 'object', this, 'fieldname', 'vertex_pairing', 'format', 'DoubleMat', 'mattype', 3);
		WriteData(fid, prefix, 'object', this, 'fieldname', 'penalty_factor', 'format', 'Double');

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
		WriteData(fid, prefix, 'data', outputs, 'name', 'md.masstransport.requested_outputs', 'format', 'StringArray');
	} //}}}

	extrude(md) {//{{{
	} //}}}
}
