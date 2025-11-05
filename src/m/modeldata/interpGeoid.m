function [geoid] = interpGeoid(X,Y,varargin),

switch oshostname(),
	case {'ronne'}
		rootname='/home/ModelData/Global/Geoid/eigen-6c4-1970.mat';
	case {'totten'}
		rootname='/totten_1/ModelData/Global/Geoid/eigen-6c4-1970.mat';
	otherwise
		error('machine not supported yet');
end
verbose = 1;

if nargin==3,
	hemisphere = varargin{1};
else
	hemisphere = +1;
end

if hemisphere==+1,
	if verbose, disp('   -- Geoid: convert to lat/lon using Greenland projection'); end
	[LAT,  LON  ] = xy2ll(double(X(:)),double(Y(:)),+1,45,70);
else
	if verbose, disp('   -- Geoid: convert to lat/lon using Antarctica projection'); end
	[LAT,  LON  ] = xy2ll(double(X(:)),double(Y(:)),-1,0,71);
end
pos=find(LON<0);
LON(pos) =360+LON(pos);

if verbose, disp('   -- Geoid: loading eigen-6c4 '); end
A=load(rootname);

if verbose, disp('   -- Geoid: interpolating'); end
geoid = InterpFromGrid(A.lon,A.lat,A.geoid,LON,LAT);
geoid = reshape(geoid,size(X));
