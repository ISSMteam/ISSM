function out = interpMartos2017(X,Y)
%INTERPMARTOS2017 - interpolate geothermal heat flux
%
%   Usage:
%      out = interpMartos2017(X,Y)

switch oshostname(),
	case {'ronne'}
		gtfpath='/home/ModelData/Antarctica/GeothermalMartos/Antarctic_GHF.xyz';
	otherwise
		error('machine not supported yet');
end

%Load data
data = load(gtfpath);

%Interpolate using nearest neighbor (dataset stops at ocean boundary!)
out = Kriging(data(:,1),data(:,2),data(:,3),X,Y,'output','nearestneighbor')/1e3; %from mW/m2 to W/m2
