function project3d() {
	//PROJECT3D - vertically project a vector from 2d mesh
	//
	//   vertically project a vector from 2d mesh (split in noncoll and coll areas) into a 3d mesh.
	//   This vector can be a node vector of size (md.mesh.numberofvertices2d,N/A) or an 
	//   element vector of size (md.mesh.numberofelements2d,N/A). 
	//   arguments: 
	//      'vector': 2d vector
	//      'type': 'element' or 'node'. 
	//   options: 
	//      'layer' a layer number where vector should keep its values. If not specified, all layers adopt the 
	//             value of the 2d vector.
	//      'padding': default to 0 (value adopted by other 3d layers not being projected0
	//
	//   Egs:
	// md.extruded_vector=project3d(md,'vector',vector2d,'type','node','layer',1,'padding',null);
	// md.extruded_vector=project3d(md,'vector',vector2d,'type','element','padding',0);
	// md.extruded_vector=project3d(md,'vector',vector2d,'type','node');

	//some regular checks
	if (arguments.length===1 || arguments.length===0) {
		console.error('project3d bad usage');
	}
	if (md.mesh.elementtype() !== 'Penta') {
		console.error('input model is not 3d');
	}

	//retrieve parameters from options.
	var options      = new pairoptions(Array.prototype.slice.call(arguments, 1)); // slice to remove md
	var vector2d     = options.getfieldvalue('vector');     //mandatory
	var type         = options.getfieldvalue('type');       //mandatory
	var layer        = options.getfieldvalue('layer',0);    //optional (do all layers default:)
	var paddingvalue = options.getfieldvalue('padding',0);  //0 by default

	if (Number.isNaN(vector2d) || vector2d === 0 || vector2d.length === 1) { // NaN treated as length 1 in MATLAB
		projected_vector=vector2d;
	} else if (type.toLowerCase() === 'node') {
		//Initialize 3d vector
		if (vector2d.length===md.mesh.numberofvertices2d) {
			projected_vector=NewArrayFill(md.mesh.numberofvertices,paddingvalue);
		} else if (vector2d.length===md.mesh.numberofvertices2d+1) {
			projected_vector=NewArrayFill(md.mesh.numberofvertices+1,paddingvalue);
			projected_vector[projected_vector.length-1] = vector2d[vector2d.length-1];
			vector2d.pop();
		} else {
			console.error('vector length not supported')
		}

		//Fill in
		if (layer===0) {
			for (var i = 1; i <= md.mesh.numberoflayers; ++i) {
				for (var j = (i-1)*md.mesh.numberofvertices2d, k = 0; j < i*md.mesh.numberofvertices2d; j++, k++) {
					projected_vector[j] = vector2d[k];
				}
			}
		} else {
			for (var j = (layer-1)*md.mesh.numberofvertices2d, k = 0; j < layer*md.mesh.numberofvertices2d; j++, k++) {
				projected_vector[j] = vector2d[k];
			}
		}
	} else if (type.toLowerCase() === 'element') {
		//Initialize 3d vector
		//var vector2d_size2 = Array.isArray(vector2d[0]) ? vector2d[0].length : 1; //get size of vector2d's 2nd axis
		if (vector2d.length===md.mesh.numberofelements2d) {
			projected_vector=NewArrayFill(md.mesh.numberofelements,paddingvalue);
		} else if (vector2d.length===md.mesh.numberofelements2d+1) {
			projected_vector=NewArrayFill(md.mesh.numberofelements+1,paddingvalue);
			projected_vector[projected_vector.length-1] = vector2d[vector2d.length-1];
			vector2d.pop();
		} else {
			console.error('vector length not supported')
		}

		//Fill in
		if (layer===0) {
			for (var i = 1; i <= md.mesh.numberoflayers-1; ++i) {
				for (var j = (i-1)*md.mesh.numberofelements2d, k = 0; j < i*md.mesh.numberofelements2d; j++, k++) {
					projected_vector[j] = vector2d[k];
				}
			}
		} else {
			for (var j = (layer-1)*md.mesh.numberofelements2d, k = 0; j < layer*md.mesh.numberofelements2d; j++, k++) {
				projected_vector[j] = vector2d[k];
			}
		}
	} else {
		console.error('project3d error message: unknown projection type');
	}

	return projected_vector;
};
