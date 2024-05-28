function plot_quiver(md, options, canvas, noCacheNodesOverride) { //{{{
	//PLOT_QUIVER - quiver plot with colors
	//
	//   Usage:
	//      plot_quiver(md, options, canvas)
	//
	//   See also: PLOTMODEL, PLOT_MANAGER

	//Disabling for now, since quivers are "sticky" - once turned on, they won't turn off. This is due to cachenodes, but should find better way to handle it.
	return;
	
	//if ('quiver' in canvas.nodes && noCacheNodesOverride && options.getfieldvalue('cachenodes','off') === 'on') return; 
	
	//{{{ declare variables:
	//Process data and model
	var meshresults = processmesh(md, [], options);
	var x = meshresults[0]; 
	var y = meshresults[1]; 
	var z = meshresults[2]; 
	var elements = meshresults[3]; 
	var is2d = meshresults[4]; 
	var isplanet = meshresults[5];
	if (!md.mesh.classname().startsWith('mesh3d')) {
		z = md.geometry.surface;
	}
		
	//Compute coordinates and data range:
	var xlim = options.getfieldvalue('xlim', [ArrayMin(x), ArrayMax(x)]);
	var ylim = options.getfieldvalue('ylim', [ArrayMin(y), ArrayMax(y)]);
	var zlim = options.getfieldvalue('zlim', [ArrayMin(z), ArrayMax(z)]);

	//Only displaying velocity fields for now
	var v = vesl.helpers.isEmptyOrUndefined(md.results) ?  md.initialization.vel : md.results[canvas.animation.frame].Vel;
	var vx = vesl.helpers.isEmptyOrUndefined(md.results) ? md.initialization.vx : md.results[canvas.animation.frame].Vx;
	var vy = vesl.helpers.isEmptyOrUndefined(md.results) ? md.initialization.vy : md.results[canvas.animation.frame].Vy;

	//Handle heightscale
	var vertices, scale;
	if (!md.mesh.classname().startsWith('mesh3d')) {
		vertices = [x, y, z];
		scale = [1, 1, options.getfieldvalue('heightscale', 1)];
	}
	else {
		vertices = Node.prototype.scaleVertices(md, x, y, z, elements, options.getfieldvalue('heightscale', 1), options.getfieldvalue('maskregion',{'enabled':false}), true);
		scale = [1, 1, 1];
	}
	
	//Compute gl variables:
	var edgecolor = options.getfieldvalue('edgecolor', 'black');
	var node = new Node(
		'canvas', canvas,
		'options', options,
		'name', 'quiver',
		'shaderName', 'Colored',
		'opacity', options.getfieldvalue('opacity', 1.0),
		//'center', [(xlim[0] + xlim[1]) / 2, (ylim[0] + ylim[1]) / 2, md.mesh.classname() === 'mesh3dsurface' ? (zlim[0] + zlim[1]) / 2 : zlim[0]],
		'center', [(xlim[0] + xlim[1]) / 2, (ylim[0] + ylim[1]) / 2, (zlim[0] + zlim[1]) / 2],
		'drawMode', canvas.gl.LINES,
		'diffuseColor', edgecolor,
		'lineWidth', options.getfieldvalue('linewidth', 1),
		'maskEnabled', options.getfieldvalue('innermask','off') == 'on',
		'maskHeight', options.getfieldvalue('innermaskheight', 150.0) / options.getfieldvalue('heightscale', 1),
		'maskColor', options.getfieldvalue('innermaskcolor',[0.0, 0.0, 1.0, 1.0]),
		'rotation', [-90, 0, 0],
		'scale', scale
	);
	
	
	//{{{ node plot
	if (elements[0].length==6){ //prisms
	}
	else if (elements[0].length==4){ //tetras
	}
	else{ //2D triangular elements
		//Create arow vertices, and use vx/vy to determine rotation before adding to quiver mesh.
		var verticesArrow = [vec3.fromValues(0.0, 0.0, 0.0), vec3.fromValues(1.0, 0.0, 0.0), vec3.fromValues(0.667, -0.167, 0.0), vec3.fromValues(1.0, 0.0, 0.0), vec3.fromValues(0.667, 0.166, 0.0), vec3.fromValues(1.0, 0.0, 0.0)];
		
		var newX = [];
		var newY = [];
		var newZ = [];
		var xyz = vec3.create();
		var vertex = vec3.create();
		var scaling = options.getfieldvalue('scaling', 1);
		var heightScale = options.getfieldvalue('heightscale', 1);
		var arrowScale;
		var modelMatrix = mat4.create();
		var scaleMatrix = mat4.create();
		var rotationMatrix = mat4.create();
		
		for(var i = 0, iX = 0, iY = 0, iZ = 0; i < x.length; i++){
			xyz = vec3.fromValues(x[i], y[i], z[i]);
			arrowScale = v[i] * scaling;
			scaleMatrix = mat4.create();
			mat4.scale(scaleMatrix, mat4.create(), vec3.fromValues(arrowScale, arrowScale, arrowScale));
			mat4.rotate(rotationMatrix, mat4.create(), Math.atan2(vy[i], vx[i]), [0.0, 0.0, 1.0]);
			mat4.multiply(modelMatrix, rotationMatrix, scaleMatrix);
			for (var j = 0; j < 6; j++){
				vec3.transformMat4(vertex, verticesArrow[j], modelMatrix);
				vec3.add(vertex, vertex, xyz);
				newX[iX++] = vertex[0];
				newY[iY++] = vertex[1];
				newZ[iZ++] = vertex[2];
			}
		}
		node.patch('Vertices', [newX, newY, newZ], 'FaceColor', 'none');
	}
	//}}}
} //}}}
