function out = interpDhdt(X,Y),

switch oshostname(),
	case {'ronne'}
		dhdtpath='/home/ModelData/Greenland/DHDT/dhdt0306.tif';
	case {'totten'}
		dhdtpath='/totten_1/ModelData/Greenland/DHDT/dhdt0306.tif';
	otherwise
		error('machine not supported yet');
end

%convert coordinates:
[lat lon] = xy2ll(X,Y,+1);
[X Y] = ll2utm(lat,lon,24);

%Get image info
Tinfo = imfinfo(dhdtpath);
N     = Tinfo.Width;
M     = Tinfo.Height;
dx    = Tinfo.ModelPixelScaleTag(1);
dy    = Tinfo.ModelPixelScaleTag(2);
minx  = Tinfo.ModelTiepointTag(4);
maxy  = Tinfo.ModelTiepointTag(5);

%Generate vectors
xdata = minx + dx/2 + ((0:N-1).*dx);
ydata = maxy - dy/2 - ((M  -1:-1:0).*dy);
ydata = fliplr(ydata);

%Get pixels we are interested in
offset=2;
xmin=min(X(:)); xmax=max(X(:));
posx=find(xdata<=xmax);
id1x=max(1,find(xdata>=xmin,1)-offset);
id2x=min(numel(xdata),posx(end)+offset);

if 0,
	ymin=min(Y(:)); ymax=max(Y(:));
	posy=find(ydata<=ymax);
	id1y=max(1,find(ydata>=ymin,1)-offset);
	id2y=min(numel(ydata),posy(end)+offset);
else
	ymin=min(Y(:)); ymax=max(Y(:));
	posy=find(ydata>=ymin);
	id1y=max(1,find(ydata<=ymax,1)-offset);
	id2y=min(numel(ydata),posy(end)+offset);
end

data  = double(imread(dhdtpath,'PixelRegion',{[id1y,id2y],[id1x,id2x]}));
xdata=xdata(id1x:id2x);
ydata=ydata(id1y:id2y);
data(find(data>+10^3)) = 0;
data(find(data<-10^3)) = 0;

out = InterpFromGrid(xdata,ydata,data,X,Y);
