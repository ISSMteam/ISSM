function [dataout] = interpBedmap(X,Y,string),
%INTERPBEDMAP - interpolate bedmap data
%
%   Available data:
%      1. bed                          is bed height
%      2. thickness                    is ice thickness
%
%   Usage:
%      [dataout] = interpBedmap(X,Y,string)

path=[jplsvn() '/proj-morlighem/DatasetAntarctica/Data/BedMap/gridded/'];

if strcmp(string,'bed'),
	path = [path '/bed.mat'];
	load(path);
	x_m =(x_m(2:end)+x_m(1:end-1))/2.;
	y_m =(y_m(2:end)+y_m(1:end-1))/2.;
	dataout = InterpFromGrid(x_m,y_m,bed,double(X),double(Y));
elseif strcmp(string,'thickness')
	path = [path '/thickness.mat'];
	load(path);
	x_m =(x_m(2:end)+x_m(1:end-1))/2.;
	y_m =(y_m(2:end)+y_m(1:end-1))/2.;
	dataout = InterpFromGrid(x_m,y_m,thickness,double(X),double(Y));
else
	error('not supported');
end
