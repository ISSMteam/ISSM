function MeshMergedOutlineToShp(mds,shapefilename)
%MeshMergedOutlineToShp - export multiple mesh outlines to shp file
%
%   Usage:
%      MeshMergedOutlineToShp({md1,md2},'Greenland.shp');

	contours= struct([]);
    count=1;
    for i=1:length(mds)
        md = mds{i};
        for j=1:md.mesh.numberofvertices
            if md.mesh.vertexonboundary(j)
                contours(count).x = md.mesh.x(j);
                contours(count).y = md.mesh.y(j);
                contours(count).id = count;
                contours(count).Geometry = 'Point';
                count = count + 1;
            end
        end
	end
	shpwrite(contours,shapefilename);
