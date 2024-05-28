'use strict';

function plot_manager(md, options, subplotwidth, nlines, ncols, i) { //{{{
	//PLOT__MANAGER - distribute the plots, called by plotmodel
	//
	//   Usage:
	//      plot_manager(md, options, subplotwidth, i);
	//
	//   See also: PLOTMODEL, PLOT_UNIT
	let canvas      = document.getElementById(options.getfieldvalue('canvasid'));
	let data        = options.getfieldvalue('data'); // Get data to be displayed
	let dataresults = null; 
	let datatype    = 0;
	let	data2       = null;
			
	// Parse options and get a structure of options
	checkplotoptions(md, options);

	// Figure out if this is a special plot
	if (typeof data === 'string') {
		switch(data) {
			case 'boundaries':
				plot_boundaries(md, options, subplotwidth, i);
				return;
				
			case 'BC':
				plot_BC(md, options, subplotwidth, i, data);
				return;
				
			case 'edges':
				plot_edges(md, options, subplotwidth, i, data);
				return;
				
			case 'elementnumbering':
				plot_elementnumbering(md, options, subplotwidth, i);
				return;
				
			case 'highlightelements':
				plot_highlightelements(md, options, subplotwidth, i);
				return;
				
			case 'qmumean':
				plot_qmumean(md, options, nlines, ncols, i);
				return;
				
			case 'qmustddev':
				plot_qmustddev(md, options, nlines, ncols, i);
				return;
				
			case 'qmuhistnorm':
				plot_qmuhistnorm(md, options, nlines, ncols, i);
				return;
				
			case 'qmu_mass_flux_segments':
				plot_qmu_mass_flux_segments(md, options, nlines, ncols, i);
				return;
				
			case 'part_hist':
				plot_parthist(md, options, nlines, ncols, i);
				return;
				
			case 'part_hist_n':
				plot_parthistn(md, options, nlines, ncols, i);
				return;
				
			case 'part_hist_w':
				plot_parthistw(md, options, nlines, ncols, i);
				return;
				
			case 'elements_type':
				plot_elementstype(md, options, subplotwidth, i);
				return;
				
			case 'vertexnumbering':
				plot_vertexnumbering(md, options, subplotwidth, i);
				return;
				
			case 'highlightvertices':
				plot_highlightvertices(md, options, subplotwidth, i);
				return;
				
			case 'basal_drag':
				plot_basaldrag(md, options, subplotwidth, i, data);
				return;
				
			case 'basal_dragx':
				plot_basaldrag(md, options, subplotwidth, i, data);
				return;
				
			case 'basal_dragy':
				plot_basaldrag(md, options, subplotwidth, i, data);
				return;
				
			case 'driving_stress':
				plot_drivingstress(md, options, subplotwidth, i);
				return;
				
			case 'mesh':
				plot_mesh(md, options, canvas);
				return;
				
			case 'none':
				return;
				
			case 'penalties':
				plot_penalties(md, options, subplotwidth, i);
				return;
				
			case 'partition':
				plot_partition(md, options, nlines, ncols, i);
				return;
				
			case 'referential':
				plot_referential(md, options, nlines, ncols, i);
				return;
				
			case 'riftvel':
				plot_riftvel(md, options, nlines, ncols, i);
				return;
				
			case 'riftnumbering':
				plot_riftnumbering(md, options, nlines, ncols, i);
				return;
				
			case 'rifts':
				plot_rifts(md, options, nlines, ncols, i);
				return;
				
			case 'riftrelvel':
				plot_riftrelvel(md, options, nlines, ncols, i);
				return;
				
			case 'riftpenetration':
				plot_riftpenetration(md, options, nlines, ncols, i);
				return;
				
			case 'riftfraction':
				plot_riftfraction(md, options, nlines, ncols, i);
				return;
				
			case 'sarpwr':
				plot_sarpwr(md, options, subplotwidth, i);
				return;
				
			case 'time_dependant':
				plot_vstime(md, options, nlines, ncols, i);
				return;
				
			case 'icefront':
				plot_icefront(md, options, subplotwidth, i, data);
				return;
				
			case 'segments':
				plot_segments(md, options, subplotwidth, i, data);
				return;
				
			case 'quiver':
				plot_quiver(md, options, canvas);
				return;
				
			case 'strainrate_tensor':
			case 'strainrate':
			case 'strainrate_principal':
			case 'strainrate_principalaxis1':
			case 'strainrate_principalaxis2':
			case 'strainrate_principalaxis3':
			case 'stress_tensor':
			case 'stress':
			case 'stress_principal':
			case 'stress_principalaxis1':
			case 'stress_principalaxis2':
			case 'stress_principalaxis3':
			case 'deviatoricstress_tensor':
			case 'deviatoricstress':
			case 'deviatoricstress_principal':
			case 'deviatoricstress_principalaxis1':
			case 'deviatoricstress_principalaxis2':
			case 'deviatoricstress_principalaxis3':
				plot_tensor(md, options, subplotwidth, i, data);
				return;
				
			case 'thermaltransient_results':
				plot_thermaltransient_results(md, options, subplotwidth, i);
				return;
				
			case 'transient_movie':
				plot_transient_movie(md, options, canvas);
				return;
				
			case 'transient_results':
				plot_transient_results(md, options, subplotwidth, i);
				return;
				
			case 'transient_field':
				plot_transient_field(md, options, subplotwidth, i);
				return;
				
			default:
				if (data in md) {
					data = md[data];
				} else {
					error('plot error message: data provided not supported yet. Type plotdoc for help');
				}
		}
	}

	// Figure out if this is a semi-transparent plot.
	if (options.getfieldvalue('overlay', 'off') === 'on') {
		plot_overlay(md, data, options, canvas);
		//return;
	}

	// Figure out if this is a semi-transparent plot.
	if (options.exist('googlemaps')) {
		plot_googlemaps(md, data, options, nlines, ncols, i);
		return;
	}

	// Figure out if this is a semi-transparent plot.
	if (options.exist('gridded')) {
		plot_gridded(md, data, options, nlines, ncols, i);
		return;
	}

	// Figure out if this is a Section plot
	if (options.exist('sectionvalue')) {
		plot_section(md, data, options, nlines, ncols, i);
		return;
	}

	//Figure out if this is a Profile plot
	if (options.exist('profile')) {
		plot_profile(md, data, options, nlines, ncols, i);
		return;
	}
	
	dataresults = processdata(md, data, options);
	data2       = dataresults[0]; 
	datatype    = dataresults[1];
	
	// Plot unit
	plot_unit(md, data2, datatype, options, canvas);

	applyoptions(md, data2, options, canvas);
} //}}}
