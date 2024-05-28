class sealevelmodel {
	/**
	 * SEALEVELMODEL class definition
	 * 
	 * Usage:
	 *     slm = sealevelmodel(varargin)
	 *
	 * Example:
	 *     slm = sealevelmodel(
	 *         'icecap', md_greenland,
	 *         'icecap', md_antarctica,
	 *         'earth', md_earth
	 *     )
	 *
	 * TODO:
	 * - Finish translation
	 */
	constructor(...varargin) {//{{{
		this.icecaps = []; // list of land/ice models; name should be changed later
		this.earth = null; // model for the whole earth
		this.basins = []; // list of basins, matching icecaps, where shapefile info is held
		this.cluster = null;
		this.miscellaneous = null;
		this.settings = null;
		this.private = null;
		this.mergedcaps = null;
		this.transitions = [];
		this.eltransitions = [];
		this.planet = '';

		this.setdefaultparamters();

		if (varargin.length) {
			let options = pairoptions(varargin);

			// Recover all the icecap models
			this.icecaps = options.getfieldvalues('ice_cap', []);

			// Recover the earth model
			this.earth = options.getfieldvalue('earth', 0);

			// Set planet type
			this.planet = options.getfieldvalue('planet', 'earth');
		}
	} //}}}

	setdefaultparamters() {//{{{
		// Initialize subclasses
		this.icecaps = [];
		this.earth = [];
		this.cluster = generic();
		this.miscellaneous = miscellaneous();
		this.settings = issmsettings();
		this.private = private();
		this.transitions = [];
		this.eltransitions = [];
		this.planet = 'earth';
	} //}}}

	static checkconsistency(slm, solutiontype) {//{{{
		// Is the coupler turned on?
		for (let i = 0; slm.icecaps.length; ++i) {
			if (!slm.icecaps[0].transient.iscoupler) {
				console.log('sealevelmodel checkconsistency error: icecap model ' +  slm.icecaps[0].miscellaneous.name + ' should have the transient coupler option turned on!');
			}
		}

		if (!slm.earth.transient.iscoupler) {
			console.log('sealevelmodel checkconsistency error: earth model should have the transient coupler option turned on!');
		}

		// Check that the transition vectors have the right size
		for (let i = 0; i < slm.icecaps.length; ++i) {
			if (slm.icecaps[0].mesh.numberorvertices != slm.earth.solidearth.transitions[0].length) {
				error('sealevelmodel checkconsistency issue with size of transition vector for ice cap: ' + i + ' name: ' + slm.icecaps[0].miscellaneous.name);
			}
		}

		// Check that runfrequency is the same everywhere
		for (let i = 0; i < slm.icecaps; ++i) {
			if (slm.icecaps[i].solidearth.settings.runfrequency != slm.earth.solidearth.settings.runfrequency) {
				error('sealevelmodel checkconsistency error: icecap model ' + slm.icecaps[i].miscellaneous.name + ' should have the same run frequency as earth!');
			}
		}

		// Make sure steric_rate is the same everywhere
		for (let i = 0; i < slm.icecaps.length; ++i) {
			let md = slm.icecaps[i];
			if (!isempty(find(md.dsl.steric_rate - slm.earth.dsl.steric_rate[slm.earth.dsl.transitions[i]]))) {
				error('steric rate on ice cap ' + md.miscellaneous.name + ' is not the same as for the earth');
			}
		}

		// Make sure grd is the same everywhere
		for (let i = 0; i < slm.icecaps.length; ++i) {
			let md = slm.icecaps[i];
			if (md.solidearthsettings.isgrd != slm.earth.solidearthsettings.isgrd) {
				error('isgrd on ice cap ' + md.miscellaneous.name + ' is not the same as for the earth');
			}
		}

		// Make sure that there is no solid earth external forcing on the basins
		for (let i = 0; i < slm.icecaps.length; ++i) {
			let md = slm.icecaps[i];
			if (!isempty(md.solidearth.external)) {
				error('cannot run external forcings on an ice sheet when running a coupling earth/ice sheet model');
			}
		}

		// Make sure that we have the right grd model for computing out sealevel patterns
		for (let i = 0; i < slm.icecaps.length; ++i) {
			let md = slm.icecaps[i];
			if (md.solidearth.settings.grdmodel) {
				error('sealevelmodel checkconsistency error message: ice sheets do not run GRD module, specify solidearth.settings.grdmodel=0 on ice cap ' + i);
			}
		}
	} //}}}

	disp() {//{{{
		console.log('WARNING: sealevelmodel::disp is not yet implemented');
	} //}}}

	mergeresults() {//{{{
		let champs = fieldnames(this.icecaps[0].results.TransientSolution);
		for (let i = 0; i < this.mergedcaps.length / 2; ++i) {
			let md = this.mergedcaps[2 * i];
			let trans = this.mergedcaps[2 * i + 1];
			//let icecaps = this.icecaps[this.range[2 * i + 2]];
			for (let j = 0; j < this.icecaps[0].results.TransientSolution) {
				for (let k = 0; champs.length; ++k) {
					if (strcmpi(typeof(icecaps[0].results.TransientSolution[j][champs[k]]) == 'double')) {
						// Vertex or element?
						if (icecaps[0].results.TransientSolution[j][champs[k]].length == icecaps[0].mesh.numberofvertices) {
							md.results.TransientSolution[j][champs[k]] == icecaps[0].mesh.numberofvertices;

						}
					} else {

					}
				}
			}
	} //}}}
}
