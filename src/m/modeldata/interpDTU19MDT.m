function mdt = interpDTU19MDT(X,Y,varargin);

switch oshostname(),
	case {'ronne'}
		rootname='/ronne_2/home/ModelData/Global/DTU19MDT/dtu19mdt.mat';
	case {'totten'}
		rootname='/totten_1/ModelData/Global/DTU19MDT/dtu19mdt.mat';
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
	if verbose, disp('   -- DTU19MDT: convert to lat/lon using Greenland projection'); end
	[LAT,  LON  ] = xy2ll(double(X(:)),double(Y(:)),+1,45,70);
else
	if verbose, disp('   -- DTU19MDT: convert to lat/lon using Antarctica projection'); end
	[LAT,  LON  ] = xy2ll(double(X(:)),double(Y(:)),-1,0,71);
end
pos=find(LON<0);
LON(pos) =360+LON(pos);
LAT=reshape(LAT,size(X));
LON=reshape(LON,size(X));

if verbose, disp('   -- DTU19MDT: loading DTU19MDT'); end
A=load(rootname);

if verbose, disp('   -- DTU19MDT: interpolating'); end
mdt = InterpFromGrid(A.lon_ext,A.lat_ext,A.mdt_ext,LON,LAT);
