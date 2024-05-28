function [bedout sid] = interpIBCSO2(X,Y),

%read data
switch (oshostname()),
	case {'totten'}
		ncpath='/totten_1/ModelData/Antarctica/IBCSO2/IBCSO_v2_bed.nc';
		sidpath='/totten_1/ModelData/Antarctica/IBCSO2/IBCSO_v2_TID.nc';
	otherwise
		error('hostname not supported yet');
end

disp('   -- IBCSOv2: Changing Coordinate system from 3031 to 9354');
[X Y]=CoordTransform(double(X),double(Y),'EPSG:3031','EPSG:9354');

disp('   -- IBCSOv2: loading bathymetry');
xdata = double(ncread(ncpath,'x'));
ydata = double(ncread(ncpath,'y'));
data  = double(ncread(ncpath,'z'))';
disp('   -- IBCSOv2: interpolating bed');
bedout = InterpFromGrid(xdata,ydata,data,double(X),double(Y));

if nargout==2,
	disp('   -- IBCSOv2: bathymetry sid');
	data  = ncread(sidpath,'tid')';
	disp('   -- IBCSOv2: interpolating sids');
	sid = InterpFromGrid(xdata,ydata,data,double(X),double(Y),'nearest');
end
