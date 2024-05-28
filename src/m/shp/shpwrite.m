function shpwrite(shp,filename)
%SHPWRITE - write a shape file from a contour structure
%
%   Usage:
%      shpwrite(shp,filename)
%
%   Example:
%      shpwrite(shp,'domainoutline.shp')
%
%   See also SHPREAD

contours=struct([]);
for i=1:length(shp),
	if strcmpi(shp(i).Geometry,'Point'),
		contours(i).Geometry='Point';
	elseif strcmpi(shp(i).Geometry,'Polygon'),
		contours(i).Geometry='Polygon';
	elseif strcmpi(shp(i).Geometry,'Line'),
		contours(i).Geometry='Line';
	else,
		error(['shpwrite error: geometry ' shp(i).Geometry ' not currently supported']);
	end
	contours(i).id=i;
	contours(i).X=shp(i).x;
	contours(i).Y=shp(i).y;
end

shapewrite(contours,filename);
end
