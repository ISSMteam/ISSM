function md=TwoDToThreeD(md,planet)
	%reproject model into lat,long if necessary:
	if ~strcmpi(md.mesh.proj,epsg2proj(4326)),
		[md.mesh.x,md.mesh.y]=CoordTransform(md.mesh.x,md.mesh.y,md.mesh.proj,'EPSG:4326');
	end

	%Make a 3dsurface mesh out of this: 
	R=planetradius(planet);

	%we assume x and y hold the long,lat values:
	long=md.mesh.x;
	lat=md.mesh.y;

	%Assume spherical body: 
	x = R .* cosd(lat) .* cosd(long);
	y = R .* cosd(lat) .* sind(long);
	z = R .* sind(lat);

	elements=md.mesh.elements;
	vc=md.mesh.vertexconnectivity;
	vb=md.mesh.vertexonboundary;
	md.mesh=mesh3dsurface();
	md.mesh.lat=lat;
	md.mesh.long=long;
	md.mesh.x=x;
	md.mesh.y=y;
	md.mesh.z=z;
	md.mesh.elements=elements;
	md.mesh.numberofelements=length(elements);
	md.mesh.numberofvertices=length(lat);
	md.mesh.r=R*ones(md.mesh.numberofvertices,1);
	md.mesh.vertexconnectivity=vc;
	md.mesh.vertexonboundary=vb;
