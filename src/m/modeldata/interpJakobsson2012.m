function [bedout sourceout] = interpJakobsson2012(X,Y,string),

switch oshostname(),
	case {'murdo','thwaites','astrid'}
		ncpath ='/u/astrid-r1b/morlighe/issmjpl/proj-morlighem/DatasetGreenland/Data/IBCAO/IBCAO_V3_500m_RR.grd';
	case {'ronne'}
		ncpath ='/home/ModelData/Greenland/IBCAO/IBCAO_V3_500m_RR.grd';
	otherwise
		error('machine not supported yet');
end

%Convert to IBCAO projections
disp('   -- Jakobsson2012: converting coordinates');
[LAT,  LON  ] = xy2ll(double(X(:)),double(Y(:)),+1,45,70);
[x0075,y0075] = ll2xy(LAT,LON,+1,0,75);

disp('   -- Jakobsson2012: loading bathymetry');
xdata = double(ncread(ncpath,'x'));
ydata = double(ncread(ncpath,'y'));
data  = double(ncread(ncpath,'z'))';

disp('   -- Jakobsson2012: interpolating bed');
bedout = InterpFromGrid(xdata,ydata,data,x0075,y0075);
bedout = reshape(bedout,size(X,1),size(X,2));

if nargout==2,
	ncpath ='/home/ModelData/Greenland/IBCAO/IBCAO_V3_SID_500m.grd';
	disp('   -- Jakobsson2012: loading source');
	xdata = double(ncread(ncpath,'x'));
	ydata = double(ncread(ncpath,'y'));
	data  = double(ncread(ncpath,'z'))';
	disp('   -- Jakobsson2012: interpolating source');
	sourceout = InterpFromGrid(xdata,ydata,data,x0075,y0075,'nearest');
	sourceout = reshape(sourceout,size(X,1),size(X,2));
end
