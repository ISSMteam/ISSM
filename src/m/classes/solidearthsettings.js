class solidearthsettings {//{{{
	/**
	 * SOLIDEARTHSETTINGS class definition
	 *
	 * Usage:
	 *     solidearthsettings = solidearthsettings();
	 */
	constructor() {//{{{
		this.reltol					= 0;
		this.abstol					= 0;
		this.maxiter				= 0;
		this.selfattraction			= 1;
		this.elastic				= 1;
		this.viscous				= 1;
		this.rotation				= 1;
		this.grdocean				= 1;
		this.ocean_area_scaling		= 0;
		this.runfrequency			= 1; // How many time steps will we skip before we run grd_core
		this.sealevelloading		= 1; // Will sea-level loads be computed?
		this.isgrd					= 0; // Will GRD patterns be computed?
		this.compute_bp_grd			= 0; // Will GRD patterns for bottom pressures be computed?
		this.degacc					= 0; // Degree increment for resolution of Green tables
		this.timeacc				= 0; // Time step accuracy required to compute Green tables
		this.horiz					= 0; // Compute horizontal deformation
		this.grdmodel				= 0; // GRD model (0 by defualt, 1 for (visco-)elastic, 2 for Ivins)
		this.cross_section_shape	= 0; // Cross section only used when GRD model is Ivins

		if (arguments.length == 0) {
			this.setdefaultparameters();
		} else {
			error('constructor not supported');
		}
	} //}}}

	disp() {//{{{
		console.log('WARNING: solidearthsettings::disp is not yet implemented');
	} //}}}

	setdefaultparameters() {//{{{
		// Convergence criterion: absolute, relative, and residual
		this.restol = 0.01; // 1 percent
		this.abstol = NaN;

		// Maximum of non-linear iterations
		this.maxiter = 5;

		// Computational flags
		this.selfattraction = 1;
		this.elastic = 1;
		this.viscous = 1;
		this.rotation = 1;
		this.grdocean = 1;
		this.ocean_area_scaling = 0;
		this.compute_bp_grd = 0;
		this.isgrd = 0;
		this.sealevelloading = 1;

		// Numerical discretiztion accuracy
		this.degacc = 0.01;
		this.timeacc = 1;

		// How many time steps we skip before we run solidearthsettings solver during transient
		this.runfrequency = 1;

		// Horizontal displacement? (not on by default)
		this.horiz = 0;

		// Cross section for Ivins model
		this.cross_section_shape = 1; // Square as default (see idege in GiaDeflectionCorex)

		// No GRD model by default
		this.grdmodel = 0;
	} //}}}

	checkconsistency(md, solution, analyses) {//{{{
		if (!analyses.includes('SealevelchangeAnalysis') || (solution === 'TransientSolution' && !md.transient.isslc)) {
			return md;
		}

		md = checkfield(md, 'fieldname', 'solidearth.settings.reltol', 'size', [1]);
		md = checkfield(md, 'fieldname', 'solidearth.settings.abstol', 'size', [1]);
		md = checkfield(md, 'fieldname', 'solidearth.settings.maxiter', 'size', [1], '>=', 1);
		md = checkfield(md, 'fieldname', 'solidearth.settings.runfrequency', 'size', [1], '>=', 1);
		md = checkfield(md, 'fieldname', 'solidearth.settings.degacc', 'size', [1], '>=', 1e-10);
		md = checkfield(md, 'fieldname', 'solidearth.settings.timeacc', 'size', [1], '>', 0);
		md = checkfield(md, 'fieldname', 'solidearth.settings.horiz', 'NaN', 1, 'Inf', 1, 'values', [0, 1]);
		md = checkfield(md, 'fieldname', 'solidearth.settings.grdmodel', '>=', 0, '<=', 2);
		md = checkfield(md, 'fieldname', 'solidearth.settings.cross_section_shape', 'numel', [1], 'values', [1, 2]);

		// Checks on computational flags
		if (this.elastic && !this.rigid) {
			error('solidearthsettings checkconsistency error message: need rigid on if elastic flag is set');
		}
		if (this.viscous && !this.elastic) {
			error('solidearthsettings checkconsistency error message: need elastic on if viscous flag is set');
		}

		// A GRD computation has been requested, make some checks on the nature of the meshes provided
		if (this.isgrd) {
			if (strcmpi(typeof(md.mesh), 'mesh3dsurface')) {
				if (this.grdmodel == 2) {
					error('model requires a 2D mesh to run gia Ivins computations (change mesh from mesh3dsurface to mesh2d)');
				}
			} else {
				if (this.grdmodel == 1) {
					error('model requires a 3D surface mesh to run GRD computations (change mesh from mesh2d to mesh3dsurface)');
				}
			}
			
			if (this.sealevelloading && !this.grdocean) {
				error('solidearthsettings checkconsistency error message: need grdocean on if sealevelloading flag is set');
			}
		}

		if (this.compute_bp_grd && !md.solidearth.settings.isgrd) {
			error('solidearthsettings checkconsistency error message; if bottom pressure grd patterns are requested, solidearth settings class should have isgrd flag on');
		}

		return md;
	} //}}}

	marshall(md, prefix, fid) {//{{{
		WriteData(fid, prefix, 'object', this, 'fieldname', 'reltol', 'name', 'md.solidearth.settings.reltol', 'format', 'Double');
		WriteData(fid, prefix, 'object', this, 'fieldname', 'abstol', 'name', 'md.solidearth.settings.abstol', 'format', 'Double');
		WriteData(fid, prefix, 'object', this, 'fieldname', 'maxiter', 'name', 'md.solidearth.settings.maxiter', 'format', 'Integer');
		WriteData(fid, prefix, 'object', this, 'fieldname', 'selfattraction', 'name', 'md.solidearth.settings.selfattraction', 'format', 'Boolean');
		WriteData(fid, prefix, 'object', this, 'fieldname', 'elastic', 'name', 'md.solidearth.settings.elastic', 'format', 'Boolean');
		WriteData(fid, prefix, 'object', this, 'fieldname', 'viscous', 'name', 'md.solidearth.settings.viscous', 'format', 'Boolean');
		WriteData(fid, prefix, 'object', this, 'fieldname', 'rotation', 'name', 'md.solidearth.settings.rotation', 'format', 'Boolean');
		WriteData(fid, prefix, 'object', this, 'fieldname', 'grdocean', 'name', 'md.solidearth.settings.grdocean', 'format', 'Boolean');
		WriteData(fid, prefix, 'object', this, 'fieldname', 'ocean_area_scaling', 'name', 'md.solidearth.settings.ocean_area_scaling', 'format', 'Boolean');
		WriteData(fid, prefix, 'object', this, 'fieldname', 'runfrequency', 'name', 'md.solidearth.settings.runfrequency', 'format', 'Integer');
		WriteData(fid, prefix, 'object', this, 'fieldname', 'degacc', 'name', 'md.solidearth.settings.degacc', 'format', 'Double');
		WriteData(fid, prefix, 'object', this, 'fieldname', 'timeacc', 'name', 'md.solidearth.settings.timeacc', 'format', 'Double', 'scale', md.constants.yts);
		WriteData(fid, prefix, 'object', this, 'fieldname', 'horiz', 'name', 'md.solidearth.settings.horiz', 'format', 'Integer');
		WriteData(fid, prefix, 'object', this, 'fieldname', 'sealevelloading', 'name', 'md.solidearth.settings.sealevelloading', 'format', 'Integer');
		WriteData(fid, prefix, 'object', this, 'fieldname', 'isgrd', 'name', 'md.solidearth.settings.isgrd', 'format', 'Integer');
		WriteData(fid, prefix, 'object', this, 'fieldname', 'compute_bp_grd', 'name', 'md.solidearth.settings.compute_bp_grd', 'format', 'Integer');
		WriteData(fid, prefix, 'object', this, 'fieldname', 'grdmodel', 'name', 'md.solidearth.settings.grdmodel', 'format', 'Integer');
		WriteData(fid, prefix, 'object', this, 'fieldname', 'cross_section_shape', 'name', 'md.solidearth.settings.cross_section_shape', 'format', 'Integer');
	} //}}}

	extrude(md) {//{{{
	} //}}}
} //}}}
