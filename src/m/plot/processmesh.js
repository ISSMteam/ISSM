function processmesh(md,data,options){ //{{{
//PROCESSMESH - process mesh to be plotted
//
//   Usage:
//      var meshresults=processmesh(md,data,options)
//      var x=meshresults[0]; 
//      var y=meshresults[1]; 
//      var z=meshresults[2]; 
//      var elements=meshresults[3]; 
//      var is2d=meshresults[4]; 
//      var isplanet=meshresults[5]; 
//
//   See also: PLOTMODEL, PROCESSDATA

	var x,y,z,elements,is2d,isplanet;

	if (md.mesh.numberofvertices==0){
		throw Error('plot error message: mesh is empty');
	}


	if (md.mesh.numberofvertices==md.mesh.numberofelements){
		throw Error(['plot error message: the number of elements is the same as the number of nodes...']);
	}

	if (options.getfieldvalue('coord','xy') !== 'latlon'){
		x=md.mesh.x.slice();
		if ('x2d' in md.mesh) x2d=md.mesh.x2d.slice();
		y=md.mesh.y.slice();
		if ('y2d' in md.mesh) y2d=md.mesh.y2d.slice();
	}
	else{
		x=md.mesh.long.slice();
		y=md.mesh.lat.slice();
	}

	if ('z' in md.mesh){
		z=md.mesh.z.slice();
	}
	else{
		z=NewArrayFill(x.length,0);
	}
	z=options.getfieldvalue('z',z);
	if (typeof z === 'string'){
		z=md[z];
	}
	
	//TODO: Make deep copy of elements array to prevent unwanted modification of model (slice creates deep copies for primitive types, shallow copies for obejcts)
	if ('elements2d' in md.mesh) elements2d=md.mesh.elements2d.slice();
	elements=md.mesh.elements.slice();

	//is it a 2d plot?
	if (md.mesh.dimension()==2){
		is2d=1;
	}
	else{
		if (options.getfieldvalue('layer',0)>=1){
			is2d=1;
		}
		else{
			is2d=0;
		}
	}

	//layer projection? 
	if (options.getfieldvalue('layer',0)>=1){
		if (options.getfieldvalue('coord','xy') === 'latlon'){
			throw Error('processmesh error message: cannot work with 3D meshes for now');
		}
		
		//we modify the mesh temporarily to a 2d mesh from which the 3d mesh was extruded. 
		x=x2d;
		y=y2d;
		z=NewArrayFill(x2d.length,0);
		elements=elements2d;
	}

	//units
	if (options.exist('unit')){
		unit=options.getfieldvalue('unit');
		x=x*unit;
		y=y*unit;
		z=z*unit;
	}

	//for now, always isplanet = 0, as we don't have the isa capability: 
	//if isa(md,'planet'),
	//	isplanet=1;
	//else
	isplanet=0;
	//end

	return  [x,y,z,elements,is2d,isplanet];
} //}}}
