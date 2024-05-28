'use strict';

/**
 * Global Function: plot_transient_movie
 *
 * Parameters:
 *		options			- Object containing a number of property/value pairs.
 *		animation		- Object containing a number of property/value pairs related to transient solution animation.
 *		progressText	- Object containing a number of property/value pairs related to updating the playback slider text on frame change.
 *		unit 			- String (default: '') which is appended to the frame index with a single white space before updating #playback-progress-text on frame change.
 *		update			- Boolean (default: true) indicating whether or not to update #playback-progress-text on frame change with the frame index. Set this to 'false' if you would rather manage the update of #playback-progress-text via 'onMovieUpdate' handler (or by some other means).
 * 
 * See Also:
 *		- /web/js/plotmodel.js
 *		- /web/js/plot_manager.js
 *
 * TODO:
 *		- Update function signature/handling to use ES6 default parameters for 'unit' and 'handling'.
 */
function plot_transient_movie(md, options, canvas) {
	//loop over the time steps
	let data = options.getfieldvalue('transient_field_data');
	let steps = new Array(data.length);
	
	for (let i = 0; i < steps.length; i++) {
		steps[i] = i;
	}
	
	//calculate caxis
	if (!options.exist('caxis')) {
		let range = [Infinity, -Infinity];
		let dataresults;
		for (let i in steps) {
			dataresults = processdata(md, data[i], options);
			range[0] = Math.min(range[0], ArrayMin(dataresults[1]));
			range[1] = Math.max(range[1], ArrayMax(dataresults[1]));
		}
		options.addfielddefault('caxis', range);
	}
	
	//Create unit node 
	let dataresults = processdata(md, data[0], options);
	let data2 = dataresults[0]; 
	let datatype = dataresults[1];
	plot_unit(md, data2, datatype, options, canvas);
	let unitOptions = options.getfieldvalue('unitOptions', {'name' : 'unit'});
	let unitName 	= unitOptions.name;
	let node 		= canvas.unitNodes[unitName];
	
	//process data
	let processedData = [];
	for (let i in steps) {
		dataresults = processdata(md, data[i].slice(), options);
		processedData[i] = dataresults[0];
	}
	
	//process mesh
	let meshresults = processmesh(md, data2, options);
	meshresults = scaleMesh(md, meshresults, options);
	let x = meshresults[0]; 
	let y = meshresults[1]; 
	let z = meshresults[2]; 
	let elements = meshresults[3];
	let is2d = meshresults[4]; 
	let isplanet = meshresults[5];
	let vertices = meshresults[6];
	let scale = meshresults[7];
	
	//process options
	if (canvas.animation !== undefined && canvas.animation.handler !== undefined) {
		clearInterval(canvas.animation.handler);
	}
	let animation = options.getfieldvalue('animation', {});
	
	let progressText = null;
	if (!vesl.helpers.isEmptyOrUndefined(animation.progressText)) {
		progressText = {
					update	: defaultFor(animation.progressText.update, true),
					unit	: defaultFor(animation.progressText.unit,  	'')
				};
	} else {
		progressText = {
					update	: false,
					unit	: ''
				};
	}
	canvas.animation = {
		frame: defaultFor(animation.frame, 			0),
		play: defaultFor(animation.play, 			true),
		fps: defaultFor(animation.fps, 				4),
		interval: defaultFor(animation.interval, 	1000 / defaultFor(animation.fps, 4)),
		loop: defaultFor(animation.loop, 			true),
		handler: 						{}, // Initialized below
		progressText: 						progressText
	}
	//display movie
	canvas.unitMovieData = processedData;
	canvas.animation.totalFrames = processedData.length;
	canvas.animation.frame = 0;

	let slider = null;
	
	if (!vesl.helpers.isEmptyOrUndefined(canvas.playbackControls)) {
		if (!vesl.helpers.isEmptyOrUndefined(canvas.playbackControls.slider)) {
			slider = canvas.playbackControls.slider;
		}
		
		
		if (!vesl.helpers.isEmptyOrUndefined(canvas.playbackControls.progressText)) {
			progressText = canvas.playbackControls.progressText;
		}
	}
	
	clearInterval(canvas.animation.handler);
	canvas.animation.handler = setInterval(function() {
		// Update current animation frame
		let frame = canvas.animation.frame;
		
		if (canvas.animation.play) {
			if (frame >= steps.length - 1) {
				if (canvas.animation.loop) {
					frame = 0;
				} else {
					canvas.animation.play = !canvas.animation.play;
				}
			} else {
				frame = (frame + 1) % steps.length;
			}
		}
		
		//If frame has changed, update unit node and data marker display.
		if (frame !== canvas.animation.lastFrame) {
			canvas.state.dispatchEvent('onMovieUpdate', canvas, frame, node);
			
			//Rescale in case surface has changed
			if (md.mesh.classname() === 'mesh3dsurface' || md.mesh.classname() === 'mesh3dprisms') {
				vertices = Node.prototype.scaleVertices(md, x, y, z, elements, options.getfieldvalue('heightscale', 1), options.getfieldvalue('maskregion',{'enabled':false, 'reproject':false}), true);
			} else {
				vertices = Node.prototype.scaleVertices(md, x, y, z, elements, options.getfieldvalue('heightscale', 1), options.getfieldvalue('maskregion',{'enabled':false, 'reproject':false}), true);
				//vertices = [x, y, z];
			}
	
			node.updateBuffer('Vertices', vertices);
			node.updateBuffer('Coords', processedData[frame]);
			canvas.unitData[unitName] = processedData[frame];
			if (canvas.graph.enabled) {
				vesl.graph.draw(canvas);
			}
			if (slider != null) {
				slider.val(frame);
			}
			if (!vesl.helpers.isEmptyOrUndefined(progressText) && canvas.animation.progressText.update) {
				if (canvas.animation.progressText.unit === '') {
					progressText.html(steps[frame].toFixed(0));
				} else {
					progressText.html(steps[frame].toFixed(0) + ' ' + canvas.animation.progressText.unit);
				}
			}
			//if (!vesl.helpers.isEmptyOrUndefined(canvas.nodes.quiver)) {
			//	plot_quiver(md,options,canvas,false);
			//}
			
		}
		
		//Save new frame info.
		canvas.animation.lastFrame = canvas.animation.frame;
		canvas.animation.frame = frame;
	}, canvas.animation.interval);
	
	//Update progress bar with new frame info.
	if (!vesl.helpers.isEmptyOrUndefined(slider)) {
		slider.max(steps.length - 1);
		slider.val(canvas.animation.frame);
	}
				
	applyoptions(md, [], options, canvas);
}
