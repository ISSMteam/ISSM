function out = interpStal2020(X,Y)
%INTERPSTAL2020 - interpolate geothermal heat flux
%
%   Usage:
%      out = interpStal2020(X,Y)

switch oshostname(),
	case {'amundsen.thayer.dartmouth.edu'}
		gtfpath='/local/ModelData/GeothermalFluxAntarcticaStal/aq1_01_20.nc';
	otherwise
		error('machine not supported yet');
end

%Load data
data = double(ncread(gtfpath,'Q'));
xcoord = double(ncread(gtfpath,'X'));
ycoord = double(ncread(gtfpath,'Y'));

%Put zero in areas with no values to avoid NaN in model
pos=find(isnan(data));
data(pos)=0;


%Interpolate using default bilinear interpolation since very coarse resolution
out = InterpFromGrid(xcoord,ycoord,data',double(X),double(Y)); %directly in W/m^2
