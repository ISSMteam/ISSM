function [output] = interpRACMO1km(X,Y),

switch oshostname(),
	case {'ronne'}
		rootname='/home/ModelData/Greenland/RACMO2_1km/SMB_MEAN1960-1989_150m.nc';
	case {'totten'}
		rootname='/totten_1/ModelData/Greenland/RACMO2_1km/SMB_MEAN1960-1989_150m.nc';
	otherwise
		error('machine not supported yet');
end
verbose = 1;

xdata = double(ncread(rootname,'xaxis'));
ydata = double(ncread(rootname,'yaxis'));

offset=2;

xmin=min(X(:)); xmax=max(X(:));
posx=find(xdata<=xmax);
id1x=max(1,find(xdata>=xmin,1)-offset);
id2x=min(numel(xdata),posx(end)+offset);

ymin=min(Y(:)); ymax=max(Y(:));
posy=find(ydata<=ymax);
id1y=max(1,find(ydata>=ymin,1)-offset);
id2y=min(numel(ydata),posy(end)+offset);

if verbose, disp('   -- RACMO 1-km: reading smb'); end
data  = double(ncread(rootname,'SMB',[id1x id1y],[id2x-id1x+1 id2y-id1y+1],[1 1]))';
xdata=xdata(id1x:id2x);
ydata=ydata(id1y:id2y);
data(find(data==-9999))=NaN;

if verbose, disp('   -- RACMO 1-km: interpolating (assuming rho_ice = 917 kg/m^3)'); end
%converting from mm / yr water eq to m/yr ice eq
data = data/1000 * 1000/917;
output = InterpFromGrid(xdata,ydata,data,double(X),double(Y));
