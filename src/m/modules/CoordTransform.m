function [xout,yout] = CoordTransform(xin,yin,projin,projout)
%COORDTRANSFORM - use PROJ to transform coordinates
%
%   Usage:
%      [xout,yout] = CoordTransform(xin,yin,projin,projout)
%      - xin,  yin : input coordinates
%      - xout, yout: output coordinates (in new projection)
%      - projin,projout: input/output projection string (PROJ)
%   
%   Examples:
%      [md.mesh.lat md.mesh.long] = CoordTransform(md.mesh.x,  md.mesh.y,   'EPSG:3413','EPSG:4326');
%      [md.mesh.x   md.mesh.y]    = CoordTransform(md.mesh.lat,md.mesh.long,'EPSG:4326','EPSG:3413');
%
%   Example of Projections:
%      lat/lon    = 'EPSG:4326'  or  lat/lon = '+proj=longlat +datum=WGS84'
%      Greenland  = 'EPSG:3413' (polar stereographic 70N 45W)
%      Antarctica = 'EPSG:3031' (polar stereographic 71S 0E)
%      IBCAO      = 'EPSG:3996' (polar stereographic 75N 0E)
%      
%   ll2xy previous default equivalent (uses with Hugues Ellispoid S)
%      Greenland  = '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.448564109 +units=m +no_defs'
%      Antarctica = '+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.448564109 +units=m +no_defs'
%      Bamber Greenland = '+proj=stere +lat_0=90 +lat_ts=71 +lon_0=-39 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
%
%   This function will only work if PROJ has been installed and --with-proj-dir
%   option has been set to its location in ISSM configuration

% Check usage
if nargin~=4
	help CoordTransform
	error('Wrong usage (see above)');
end

% If this function is called from within one of our distributable packages, set 
% the 'PROJ_LIB' environment variable so that the PROJ binary can find the 
% assets it needs
whatproj=what('share/proj');
if ~isempty(whatproj) && isdir(whatproj.path)
	setenv('PROJ_LIB', whatproj.path);
elseif exist([issmdir() 'externalpackages/proj/install/share/proj/'], 'dir')
	setenv('PROJ_LIB', [issmdir() 'externalpackages/proj/install/share/proj/']);
end
setenv('DYLD_LIBRARY_PATH',[issmdir() '/externalpackages/proj/install/lib']);

% Call mex module
[xout, yout] = CoordTransform_matlab(xin,yin,projin,projout);

