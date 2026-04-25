function [vxout vyout] = interpJoughinCompositeGreenland(X,Y,path)
% INTERPJOUGHINCOMPOSITEGREENLAND - interpolate Joughin's mosaic nsidc-0670
%
%   Usage:
%      [vx vy] = interpJoughinCompositeGreenland(X,Y)
%      [vx vy] = interpJoughinCompositeGreenland(X,Y,path)
%      vel     = interpJoughinCompositeGreenland(X,Y)
%
%   Example:
%      [vx vy] = interpJoughinCompositeGreenland(md.mesh.x, md.mesh.y)
%      [vx vy] = interpJoughinCompositeGreenland(md.mesh.x, md.mesh.y, '../Data')

%possible paths of dataset
paths = {...
	['/totten_1/ModelData/Greenland/VelMEaSUREs/Greenland_1995_2015_decadal_average_mosaic_v1/',],...
	['/home/ModelData/Antarctica/VelMEaSUREs/Greenland_1995_2015_decadal_average_mosaic_v1/',],...
	[issmdir() 'examples/Data/'],...
	['./',],...
	};
if nargin>2
	paths{end+1} = path;
end

%Check if we can find it
found = 0;
for i=1:numel(paths)
	if exist([paths{i} '/greenland_vel_mosaic250_vx_v1.tif'],'file')
		datadir = paths{i};
		found = 1;
		break;
	end
end
if ~found
	error(['Could not find greenland_vel_mosaic250_vx_v1.tif, you can add the path to the list or provide its path as a 3rd argument']);
end

%Now go ahead and interpolate
vxout = interpFromGeotiff([datadir '/greenland_vel_mosaic250_vx_v1.tif'], X, Y,-2e9);
vyout = interpFromGeotiff([datadir '/greenland_vel_mosaic250_vy_v1.tif'], X, Y,-2e9);

%return vel if nargount ==1
if nargout==1
	vxout = sqrt(vxout.^2+vyout.^2);
end
