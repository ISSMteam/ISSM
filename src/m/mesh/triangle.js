function triangle(md){
//TRIANGLE - create model mesh using the triangle package
//
//   This routine creates a model mesh using Triangle and a domain outline, to within a certain resolution
//   where md is a @model object, domainname is the name of an Argus domain outline file, 
//   and resolution is a characteristic length for the mesh (same unit as the domain outline
//   unit). Riftname is an optional argument (Argus domain outline) describing rifts.
//
//   Usage:
//      triangle(md,domain,resolution)
//   or triangle(md,domain,riftname, resolution)
//
//   Examples:
//      triangle(md,domain,1000);
//      triangle(md,domain, rifts, 1500);

	if (!(arguments.length==3 | arguments.length==4)){
		console.log('triangle usage error.');
	}
	
	var md=arguments[0];
	var domain=arguments[1];

	if (arguments.length==3){
		var resolution=arguments[2];
		var rifts=[];
	}
	if (arguments.length==4){
		var rifts=arguments[2];
		var resolution=arguments[3];
	}

	//Figure out a characteristic area. Resolution is a node oriented concept (ex a 1000m  resolution node would 
	//be made of 1000*1000 area squares). 
	var area=Math.pow(resolution,2);

	//Call mesher: 
	var return_array=Triangle(md, domain, rifts, area); 

	//Plug into md:
	md.mesh.elements=return_array[0];
	md.mesh.x=return_array[1];
	md.mesh.y=return_array[2];
	md.mesh.segments=return_array[3];
	md.mesh.segmentmarkers=return_array[4];
	
	//Fill in rest of fields:
	md.mesh.numberofelements=md.mesh.elements.length;
	md.mesh.numberofvertices=md.mesh.x.length;
	md.mesh.vertexonboundary=new Float64Array(md.mesh.numberofvertices); 

	for (i=0;i<md.mesh.segments.length;i++) for(var j=0;j<2;j++) md.mesh.vertexonboundary[md.mesh.segments[i][j]-1]=1;

	//Now, build the connectivity tables for this mesh.
	md.mesh.vertexconnectivity=NodeConnectivity(md.mesh.elements,md.mesh.numberofvertices);
	md.mesh.elementconnectivity=ElementConnectivity(md.mesh.elements,md.mesh.vertexconnectivity);	

}
