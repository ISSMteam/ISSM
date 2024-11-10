function mdt = interpDTU19MDT(X,Y, hemisphere);

switch oshostname(),
	case {'ronne'}
		rootname='/ronne_2/home/ModelData/Global/DTU19MDT/dtu19mdt.mat';
	case {'totten'}
		rootname='/totten_1/ModelData/Global/DTU19MDT/dtu19mdt.mat';
	otherwise
		error('machine not supported yet');
end
verbose = 1;


if hemisphere==+1,
	if verbose, disp('   -- DTU19MDT: convert to lat/lon using Greenland projection'); end
	[LAT,  LON  ] = xy2ll(double(X(:)),double(Y(:)),+1,45,70);
elseif hemisphere==-1;
	if verbose, disp('   -- DTU19MDT: convert to lat/lon using Antarctica projection'); end
	[LAT,  LON  ] = xy2ll(double(X(:)),double(Y(:)),-1,0,71);
else
	error('not supported yet');
end
pos=find(LON<0);
LON(pos) =360+LON(pos);
LAT=reshape(LAT,size(X));
LON=reshape(LON,size(X));

if verbose, disp('   -- DTU19MDT: loading DTU19MDT'); end
A=load(rootname);

if verbose, disp('   -- DTU19MDT: interpolating'); end
mdt = InterpFromGrid(A.lon_ext,A.lat_ext,A.mdt_ext,LON,LAT);
