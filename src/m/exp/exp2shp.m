function exp2shp(expfilename,shpfilename,geometry)
%SHPWRITE - write a shape file from a contour structure
%
%   Usage:
%      exp2shp(expfilename,shpfilename,geometry)
%
%   Example:
%      exp2shp('domainoutline.exp','domainoutline.shp')
%      exp2shp('domainoutline.exp','domainoutline.shp','Polygon')
%      exp2shp('massfluxgate.exp','massfluxgate.shp','Line')
%
%   See also SHPREAD, SHPWRITE, SHP2EXP

%check file extensions
[pathstr,name,ext] = fileparts(shpfilename);
if ~strcmp(ext,'.shp'),
	error(['Shapefile ' shpfilename ' does not have an extension .shp']);
end

[pathstr,name,ext] = fileparts(expfilename);
if ~strcmp(ext,'.exp'),
	error(['Exp file ' expfilename ' does not have an extension .exp']);
end

shp=expread(expfilename);

%initialize number of profile
count=1;

contours=struct([]);
for i=1:length(shp),
	if nargin < 3

		%TEMP
		%if contains(shp(i).name,'_pointcloud');
		%	continue;
		%end

		if length(shp(i).x) == 0
			continue;
		elseif contains(shp(i).name,'_pointcloud');
			geometry = 'MultiPoint';
			shp(i).name = erase(shp(i).name,'_pointcloud');
		elseif length(shp(i).x) == 1
			geometry = 'Point';
		elseif length(shp(i).x) < 3
			geometry = 'Line';
		else 
			if (shp(i).x(end)==shp(i).x(1) && shp(i).y(end)==shp(i).y(1)),
				geometry = 'Polygon';
			else
				geometry = 'Line';
			end
		end
	end
	contours(count).Geometry=geometry;
	contours(count).id=i;
	contours(count).Name=shp(i).name;
	contours(count).X=shp(i).x;
	contours(count).Y=shp(i).y;
	count = count+1;
end

%Make sure it is one single geometry otherwise it will yell at you
shapewrite(contours,shpfilename);
