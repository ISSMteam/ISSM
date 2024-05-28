'use strict';

function plotmodel(md) { //{{{
	let args            = Array.prototype.slice.call(arguments); // Convert arguments to array 
    let canvas          = null;
    let ncols           = 0;
    let nlines          = 0;
    let numberofplots   = 0;
    let options         = new plotoptions(args.slice(1, args.length)); // Process options
    let subplotwidth    = Math.ceil(Math.sqrt(options.numberofplots)); // Get number of subplots
    
	// Get figure number and number of plots
	numberofplots = options.numberofplots;

	// If nlines and ncols specified, then bypass.
	if (options.list[0].exist('nlines')) {
		nlines = options.list[0].getfieldvalue('nlines');
	} else {
		nlines = Math.ceil(numberofplots / subplotwidth);
	}
	
	if (options.list[0].exist('ncols')){
		ncols = options.list[0].getfieldvalue('ncols');
	} else {
		ncols = subplotwidth;
	}
	
	//check that nlines and ncols were given at the same time!
	if ((options.list[0].exist('ncols') && !options.list[0].exist('nlines')) 
	    || ((options.list[0].exist('nlines') && !options.list[0].exist('ncols')))) {
    	    throw Error('plotmodel error message: nlines and ncols need to be specified together, or not at all');
    }

	// Go through subplots
	if (numberofplots) {
		// Go through all data plottable and close window if an error occurs

		let canvases = [];
		for (let i = 0; i < numberofplots; ++i) {
			canvases[i] = document.getElementById(options.list[i].getfieldvalue('canvasid'));
			canvases[i].state.dispatchEvent('onPlotModelStart', canvases[i], options.list[i]);
		}
		for (let i = 0; i < numberofplots; ++i) {
			canvases[i].state.dispatchEvent('onPlotManagerStart', canvases[i], options.list[i]);
			plot_manager(options.list[i].getfieldvalue('model',md), options.list[i], subplotwidth, nlines, ncols, i);
			canvases[i].state.dispatchEvent('onPlotManagerEnd', canvases[i], options.list[i]);

			// List all unused options
			options.list[i].displayunused();
		}
		for (let i = 0; i < numberofplots; ++i) {
			canvases[i].state.dispatchEvent('onPlotModelEnd', canvases[i], options.list[i]);
		}
	}
} //}}}
