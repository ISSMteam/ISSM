function [bedout sid] = interpIBCSO(X,Y),

%read data
switch (oshostname()),
	case {'ronne'}
		ncpath='/home/ModelData/Antarctica/IBCSO/ibcso_v1_bed.grd';
		sidpath='/home/ModelData/Antarctica/IBCSO/ibcso_v1_sid.grd';
	case {'totten'}
		ncpath='/totten_1/ModelData/Antarctica/IBCSO/ibcso_v1_bed.grd';
		sidpath='/totten_1/ModelData/Antarctica/IBCSO/ibcso_v1_sid.grd';
	otherwise
		error('hostname not supported yet');
end

disp('   -- IBCSO: loading bathymetry');
x_range = double(ncread(ncpath,'x_range'));
y_range = double(ncread(ncpath,'y_range'));
spacing = double(ncread(ncpath,'spacing'));
xdata = (x_range(1)-spacing(1)/2) : spacing(1) : (x_range(2)-spacing(1)/2); 
ydata = (y_range(1)-spacing(2)/2) : spacing(2) : (y_range(2)-spacing(2)/2); 
data  = double(ncread(ncpath,'z'));
data(find(data==-9999 | isinf(data))) = NaN;
data  = reshape(data,[numel(xdata) numel(ydata)])';
disp('   -- IBCSO: interpolating bed');
bedout = InterpFromGrid(xdata,fliplr(ydata),data,double(X),double(Y));

if nargout==2,
	disp('   -- IBCSO: bathymetry sid');
	xdata = ncread(sidpath,'x');
	ydata = ncread(sidpath,'y');
	data  = ncread(sidpath,'z')';
	disp('   -- IBCSO: transforming coordinates');
	[LAT,LON] = xy2ll(double(X(:)),double(Y(:)),-1,0,71);
	[x065,y065] = ll2xy(LAT,LON,-1,0,65);
	x065 = reshape(x065,size(X));
	y065 = reshape(y065,size(Y));
	disp('   -- IBCSO: interpolating sids');
	sid = InterpFromGrid(xdata,ydata,data,x065,y065,'nearest');
	sid(find(sid<200000)) = 0;
	sid(find(sid>399999)) = 0;
end
