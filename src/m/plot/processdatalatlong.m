function [data datatype] = processdatalatlong(md,data,options);

	%   datatype = 1 -> elements
	%   datatype = 2 -> nodes
	%what is the mesh we are using: 
	x0=md.mesh.long;
	y0=md.mesh.lat;

	%add row at lat=90 and lat=-90
	add=[(-180:.1:-1e-5)';(1e-5:.1:180)'];
	nadd=length(add);
	xextra=[add;add];
	yextra=[90*ones(nadd,1); -90*ones(nadd,1)];
	x=[x0;xextra];
	y=[y0;yextra];
	elements=delaunay(x,y);
	
	%with this mesh, interpolate data: 
	if length(data)==length(md.mesh.long),
		datatype=2;

		%interpolate data: 
		extradata=griddata(x0,y0,data,xextra,yextra,'nearest');
		data=[data; extradata];
	elseif length(data)==length(md.mesh.elements),
		datatype=1;
	end
