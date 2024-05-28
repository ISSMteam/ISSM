function plot_mesh(md, options, canvas) { //{{{
	//PLOT_MESH - Function for plotting wireframe mesh.
	//
	//   Usage:
	//      plot_mesh(md, options, canvas);
	//
	//   See also: PLOTMODEL, PLOT_MANAGER

	//{{{
	//Process data and model
	let meshresults = processmesh(md, [], options);
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
	let edgecolor = options.getfieldvalue('edgecolor', 'black');
	let node = new Node(
		'canvas', canvas,
		'options', options,
		'md', md,
		'name', 'mesh',
		'cullFace', THREE.DoubleSide,
		'shaderName', 'lines',
		'opacity', options.getfieldvalue('opacity', 1.0),
		'lineWidth', options.getfieldvalue('linewidth', 1),
		'scale', scale
	);
	//}}}
	//{{{ node plot
	if (elements[0].length === 6) { //prisms
		//We can skip bottom and every other side to avoid drawing edges twice.
		let abc = elements.map(function(value, index) { return [value[0], value[1], value[2]]; }); //top
		// let dfe = elements.map(function(value, index) { return [value[3], value[5], value[4]]; }); //bottom
		// let aeb = elements.map(function(value, index) { return [value[0], value[4], value[1]]; }); //1st side upper right
		let ade = elements.map(function(value, index) { return [value[0], value[3], value[4]]; }); //1st side lower left
		// let bfc = elements.map(function(value, index) { return [value[1], value[5], value[2]]; }); //2nd side upper right
		let bef = elements.map(function(value, index) { return [value[1], value[4], value[5]]; }); //2nd side lower left
		// let cda = elements.map(function(value, index) { return [value[2], value[3], value[0]]; }); //3rd side upper right
		let cfd = elements.map(function(value, index) { return [value[2], value[5], value[3]]; }); //3rd side lower left
		let prismElements = abc.concat(ade, bef, cfd);
		node.patch('Faces', prismElements, 'Vertices', vertices, 'FaceColor', 'none', 'EdgeColor', edgecolor);
	} else if (elements[0].length === 4) { //tetras
        // TODO: Implement handling for tetras
	} else { // 2D triangular elements
		node.patch('Faces', elements, 'Vertices', vertices, 'FaceColor', 'none', 'EdgeColor', edgecolor);
	}
	//}}}
	//options=options.addfielddefault('title','Mesh');
	//options=addfielddefault('colorbar','off');
	//applyoptions(md,[],options,canvas);
} //}}}
