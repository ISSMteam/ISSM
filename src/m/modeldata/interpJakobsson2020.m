function [bedout sourceout] = interpJakobsson2020(X,Y,string),

switch oshostname(),
	case {'ronne'}
		ncpath ='/home/ModelData/Greenland/IBCAO/IBCAO_v4_200m.nc';
	case {'totten'}
		ncpath ='/totten_1/ModelData/Greenland/IBCAO/IBCAO_v4_200m.nc';
	otherwise
		error('machine not supported yet');
end

%Convert to IBCAO projections
disp('   -- Jakobsson2020: converting coordinates');
[LAT,  LON  ] = xy2ll(double(X(:)),double(Y(:)),+1,45,70);
[x0075,y0075] = ll2xy(LAT,LON,+1,0,75);

disp('   -- Jakobsson2020: loading coordinates');
xdata = double(ncread(ncpath,'x'));
ydata = double(ncread(ncpath,'y'));

offset=2;

xmin=min(x0075(:)); xmax=max(x0075(:));
posx=find(xdata<=xmax);
id1x=max(1,find(xdata>=xmin,1)-offset);
id2x=min(numel(xdata),posx(end)+offset);

ymin=min(y0075(:)); ymax=max(y0075(:));
posy=find(ydata>=ymin);
id1y=max(1,find(ydata<=ymax,1)-offset);
id2y=min(numel(ydata),posy(end)+offset);

disp(['   -- Jakobsson2020: loading bathymetry']);
data = double(ncread(ncpath,'z',[id1x id1y],[id2x-id1x+1 id2y-id1y+1],[1 1]))';
xdata=xdata(id1x:id2x);
ydata=ydata(id1y:id2y);

disp('   -- Jakobsson2020: interpolating bed');
bedout = InterpFromGrid(xdata,ydata,data,x0075,y0075);
bedout = reshape(bedout,size(X,1),size(X,2));

if nargout==2,
	ncpath ='/totten_1/ModelData/Greenland/IBCAO/IBCAO_v4_200m_TID.nc';
	disp('   -- Jakobsson2020: loading source');
	data = double(ncread(ncpath,'z',[id1x id1y],[id2x-id1x+1 id2y-id1y+1],[1 1]))';
	disp('   -- Jakobsson2020: interpolating source');
	sourceout = InterpFromGrid(xdata,ydata,data,x0075,y0075,'nearest');
	sourceout = reshape(sourceout,size(X,1),size(X,2));
end
