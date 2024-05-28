function MeshOutlineToShp(md,shapefilename)
%MESHOUTLINETOSHP - export mesh outline to shp file
%
%   Usage:
%      MeshOutlineToShp(md,'Greenland.shp');

	contours= struct([]);
	for i=1:md.mesh.numberofvertices
        if md.mesh.vertexonboundary(i)
            contours(i).x = md.mesh.x(i);
            contours(i).y = md.mesh.y(i);
            contours(i).id = i;
            contours(i).Geometry = 'Point';
        end
	end
	shpwrite(contours,shapefilename);
