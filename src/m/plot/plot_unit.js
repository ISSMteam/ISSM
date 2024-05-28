'use strict';

function plot_unit(md, data, datatype, options, canvas) { //{{{
	//PLOT_UNIT - unit plot, display data
	//
	//   Usage:
	//      plot_unit(md, data, options, canvas);
	//
	//   See also: PLOTMODEL, PLOT_MANAGER
	
	//{{{
	// Process data and model
	let meshresults = processmesh(md, data, options);
	meshresults = scaleMesh(md, meshresults, options);
	let x = meshresults[0]; 
	let y = meshresults[1]; 
	let z = meshresults[2]; 
	let elements = meshresults[3];
	let is2d = meshresults[4]; 
	let isplanet = meshresults[5];
	let vertices = meshresults[6];
	let scale = meshresults[7];
	
	//Compute gl variables:
	let edgecolor = options.getfieldvalue('edgecolor', [1.0, 1.0, 1.0, 1.0]);
	let maskzeros = options.getfieldvalue('maskzeros', {});
	
	let unitOptions = options.getfieldvalue('unitOptions', {'name' : 'unit'});
	let unitName 	= unitOptions.name;
	
	let node = new Node(
		'canvas', canvas,
		'options', options,
		'md', md,
		'name', unitName,
		'shaderName', 'unit',
		'opacity', options.getfieldvalue('opacity', 1.0),
		'caxis', options.getfieldvalue('caxis',[ArrayMin(data), ArrayMax(data)]),
		'cullFace', THREE.FrontSide,
		'enabled', options.getfieldvalue('nodata', 'off') === 'off',
		'log', options.getfieldvalue('log',false),
		'maskObject', options.getfieldvalue('maskregion', {'enabled': false}),
		'maskZerosColor', defaultFor(maskzeros.color, [1.0, 1.0, 1.0, 1.0]),
		'maskZerosEnabled', defaultFor(maskzeros.enabled, false),
		'maskZerosTolerance', defaultFor(maskzeros.tolerance, 1e-3),
		'maskZerosZeroValue', defaultFor(maskzeros.zeroValue, 0.5),
		'scale', scale
	);
	//}
	
	if (vesl.helpers.isEmptyOrUndefined(canvas.unitNodes)) {
		canvas.unitNodes = {};
	}
	
	if (vesl.helpers.isEmptyOrUndefined(canvas.unitData)) {
		canvas.unitData = {};
	}
	
	// If option 'clf' is on, remove all cached data associated with this unit
	if (options.getfieldvalue('clf', 'on') === 'on') {
		delete canvas.unitNodes[unitName];
		canvas.unitData[unitName] = data;
		
		for (let i = 0, numUnits = canvas.state.scene.children.length; i < numUnits; ++i) {
			if (canvas.state.scene.children[i].name === unitName) {
				canvas.state.scene.children.splice(i, 1);
				break;
			}
		}
	}
	
	canvas.unitNodes[unitName] = node;
	
	//}}}
	switch(datatype) {
		//{{{ element plot
		case 1:
			//WARNING: NaN are not properly found (NaN != NaN = true)
			let pos = ArrayFindNot(data, NaN); //needed for element on water
			
			if (elements[0].length === 6) { // prisms
				let abc = elements.map(function(value, index) { return [value[0], value[1], value[2]]; }); //top
				let dfe = elements.map(function(value, index) { return [value[3], value[5], value[4]]; }); //bottom
				let aeb = elements.map(function(value, index) { return [value[0], value[4], value[1]]; }); //1st side upper right
				let ade = elements.map(function(value, index) { return [value[0], value[3], value[4]]; }); //1st side lower left
				let bfc = elements.map(function(value, index) { return [value[1], value[5], value[2]]; }); //2nd side upper right
				let bef = elements.map(function(value, index) { return [value[1], value[4], value[5]]; }); //2nd side lower left
				let cda = elements.map(function(value, index) { return [value[2], value[3], value[0]]; }); //3rd side upper right
				let cfd = elements.map(function(value, index) { return [value[2], value[5], value[3]]; }); //3rd side lower left
				let prismElements = abc.concat(dfe, aeb, ade, bfc, bef, cda, cfd);
				node.patch('Faces', prismElements, 'Vertices', vertices, 'FaceVertexCData', data, 'FaceColor', 'flat', 'EdgeColor', edgecolor);
			} else if (elements[0].length === 4) { // tetras
    			// TODO: Implement handling for tetras
			} else { // triangular elements
				node.patch('Faces', elements, 'Vertices', vertices, 'FaceVertexCData', data, 'FaceColor', 'flat', 'EdgeColor', edgecolor);
			}
			
			break;
		//}}}
		//{{{ node plot
		case 2:
			if (elements[0].length === 6) { // prisms
				let abc = elements.map(function(value, index) { return [value[0], value[1], value[2]]; }); //top
				let dfe = elements.map(function(value, index) { return [value[3], value[5], value[4]]; }); //bottom
				let aeb = elements.map(function(value, index) { return [value[0], value[4], value[1]]; }); //1st side upper right
				let ade = elements.map(function(value, index) { return [value[0], value[3], value[4]]; }); //1st side lower left
				let bfc = elements.map(function(value, index) { return [value[1], value[5], value[2]]; }); //2nd side upper right
				let bef = elements.map(function(value, index) { return [value[1], value[4], value[5]]; }); //2nd side lower left
				let cda = elements.map(function(value, index) { return [value[2], value[3], value[0]]; }); //3rd side upper right
				let cfd = elements.map(function(value, index) { return [value[2], value[5], value[3]]; }); //3rd side lower left
				let prismElements = abc.concat(dfe, aeb, ade, bfc, bef, cda, cfd);
				node.patch('Faces', prismElements, 'Vertices', vertices, 'FaceVertexCData', data, 'FaceColor', 'interp', 'EdgeColor', edgecolor);
			} else if (elements[0].length === 4) { // tetras
    			// TODO: Implement handling for tetras
			} else { // triangular elements	
				node.patch('Faces', elements, 'Vertices', vertices, 'FaceVertexCData', data, 'FaceColor', 'interp', 'EdgeColor', edgecolor);
			}
			
			break;
		//}}}
		//{{{ quiver plot 
		case 3:
			if (is2d) {
				//plot_quiver(x, y, data(:, 1), data(:, 2), options);
			} else {
				//plot_quiver3(x, y, z, data(:, 1), data(:, 2), data(:, 3), options);
			}
			
			break;
		//}}}
		default:
			throw Error(sprintf('%s%i%s\n','case ', datatype,' not supported'));
	}
} //}}}
