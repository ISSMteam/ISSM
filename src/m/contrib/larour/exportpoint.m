function exportpoint(x,y,shapefilename); 

	contours=struct([]);
	for i=1:length(x),
		contours(i).id=i;
		contours(i).Lon=x(i);
		contours(i).Lat=y(i);
		contours(i).Geometry='Point';
	end
	shapewrite(contours,shapefilename);
