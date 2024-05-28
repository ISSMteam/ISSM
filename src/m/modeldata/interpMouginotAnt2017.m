function [vxout vyout]= interpMouginotAnt2017(X,Y),

%read data
switch (oshostname()),
	case {'ronne'}
		nc = '/home/ModelData/Antarctica/MouginotVel/vel_nsidc.CF16_2.nc';
	case {'totten'}
		nc = '/totten_1/ModelData/Antarctica/MouginotVel/vel_nsidc.CF16_2.nc';
	otherwise
		error('hostname not supported yet');
end

xdata = double(ncread(nc,'x'));
ydata = double(ncread(nc,'y'));

offset=2;

xmin=min(X(:)); xmax=max(X(:));
posx=find(xdata<=xmax);
id1x=max(1,find(xdata>=xmin,1)-offset);
id2x=min(numel(xdata),posx(end)+offset);

ymin=min(Y(:)); ymax=max(Y(:));
posy=find(ydata>=ymin);
id1y=max(1,find(ydata<=ymax,1)-offset);
id2y=min(numel(ydata),posy(end)+offset);

disp(['   -- Mouginot 2017: loading velocities']);
vxdata = double(ncread(nc,'VX',[id1x id1y],[id2x-id1x+1 id2y-id1y+1],[1 1]))';
vydata = double(ncread(nc,'VY',[id1x id1y],[id2x-id1x+1 id2y-id1y+1],[1 1]))';
xdata=xdata(id1x:id2x);
ydata=ydata(id1y:id2y);

disp(['   -- Mouginot 2017: interpolating ']);
vxout = InterpFromGrid(xdata,ydata,vxdata,double(X),double(Y));
vyout = InterpFromGrid(xdata,ydata,vydata,double(X),double(Y));

%return vel if only one output is requested
if nargout==1,
	vxout = sqrt(vxout.^2+vyout.^2);
end
