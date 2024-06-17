//GMSHPLANET - mesh generation for a sphere. Very specific code for gmsh from $ISSM_DIR/src/demos/simple_geo/sphere.geo
//
//   Available options (for more details see ISSM website http://issm.jpl.nasa.gov/):
//
//   - radius:             radius of the planet in km
//   - resolution:         resolution in km
//   - refine:             provide mesh
//   - refinemetric:       mesh quantity to specify resolution
//
//   Returns 'mesh3dsurface' type mesh
//
//   Examples:
//      md.mesh=gmshplanet('radius',6000,'resolution',100);
//      md.mesh=gmshplanet('radius',6000,'resolution',100);
async function gmshplanet(md) {
	let args = Array.prototype.slice.call(arguments, 1);
	let options = new pairoptions(args);
	let radius = options.getfieldvalue('radius', planetradius('earth')); //warning: planetradius in m, not km
	let resolution = options.getfieldvalue('resolution', 700);

	md.cluster = cluster;
	md = await solve(md, 'gmsh', 'checkconsistency', 'false', 'hmin', resolution);

	//methods 
	let results = md.results[0];
	let mesh = new mesh3dsurface();
	let elements = [];
	for (let i = 0; i < results.Elements.length; i=i+3){
		elements.push(results.Elements.slice(i, i+3));
	}
	mesh.numberofelements = results.NumberOfElements;
	mesh.numberofvertices =  results.NumberOfVertices;
	mesh.x = results.X;
	mesh.y = results.Y;
	mesh.z = results.Z;
	mesh.r = results.R;
	mesh.elements = elements;
	mesh.lat = results.Lat;
	mesh.long = results.Long;
	mesh.average_vertex_connectivity = results.AverageVertexConnectivity;

	//a little technicality here. the mesh generate is not exactly on the 
	//sphere. we create lat,long coordinates, and reproject on an exact sphere. 
	mesh.r=ArraySqrt(ArrayAdd(ArrayAdd(ArrayPow(mesh.x,2),ArrayPow(mesh.y,2)),ArrayPow(mesh.z,2)));

	//make sure we don't have south and north pole: 
	let pos=find(ArrayAnd(ArrayEqual(mesh.x,0), ArrayEqual(mesh.y,0)));
	mesh.lat = asind(ArrayDivide(mesh.z,mesh.r));
	mesh.long = atan2d(mesh.y,mesh.x);
	pos=find(ArrayEqual(mesh.lat,90)); ArrayIndex(mesh.lat,pos,90-.01);
	pos=find(ArrayEqual(mesh.lat,-90)); ArrayIndex(mesh.lat,pos,-90+.01);

	let radius = planetradius('earth');
	mesh.r=NewArrayFill(mesh.numberofvertices,radius);
	mesh.x=ArrayMultiply(radius,ArrayMultiply(cosd(mesh.lat),cosd(mesh.long)));
	mesh.y=ArrayMultiply(radius,ArrayMultiply(cosd(mesh.lat),sind(mesh.long)));
	mesh.z=ArrayMultiply(radius,sind(mesh.lat));

	return mesh;
}
