function MeshToShp(md,shapefilename)
%MESHTOSHP - export mesh to shp file
%
%   Usage:
%      MeshToShp(md,'Greenland.shp');

	contours= struct([]);
	for i=1:md.mesh.numberofvertices,
		contours(i).x = md.mesh.x(i);
		contours(i).y = md.mesh.y(i);
		contours(i).id = i;
		contours(i).Geometry = 'Point';
	end
	shpwrite(contours,shapefilename);
