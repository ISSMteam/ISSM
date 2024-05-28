function temp = interpISMIP6Temp(X, Y)
%interpISMIP6Temp - interpolate land ice basal temperature to the given mesh
%
%	X and Y are the coordinates of the mesh
%
filename = '/totten_1/ModelData/Greenland/ISMIP6/GreenlandISMIP6-Morlighem-2020-10-01.nc';

x = ncread(filename, 'x');
y = ncread(filename, 'y');
tb = ncread(filename, 'tb');

temp = InterpFromGrid(double(x), double(y), double(tb'), X, Y);
