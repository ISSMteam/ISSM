function sout = interpREMA(X,Y),

switch oshostname(),
	case {'ronne'}
		remapath='/home/ModelData/Antarctica/REMA/REMA_200m_dem_filled.tif';
	case {'totten'}
		remapath='/totten_1/ModelData/Antarctica/REMA/REMA_200m_dem_filled.tif';
	otherwise
		error('machine not supported yet');
end

usemap = 0;
if license('test','map_toolbox')==0,
	disp('WARNING: map toolbox not installed, trying house code');
	usemap = 0;
elseif license('checkout','map_toolbox')==0
	disp('WARNING: map toolbox not available (checkout failed), trying house code');
	usemap = 0;
end

if usemap,
	[data,R] = geotiffread(remapath);
	data=double(flipud(data));
	xdata=R.XLimWorld(1):R.DeltaX:R.XLimWorld(2); xdata=xdata(:);
	xdata =(xdata(1:end-1)+xdata(2:end))/2;
	ydata=R.YLimWorld(2):R.DeltaY:R.YLimWorld(1); ydata=flipud(ydata(:));
	ydata =(ydata(1:end-1)+ydata(2:end))/2;
else

	%Get image info
	Tinfo = imfinfo(remapath);
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

	ymin=min(Y(:)); ymax=max(Y(:));
	posy=find(ydata>=ymin);
	id1y=max(1,find(ydata<=ymax,1)-offset);
	id2y=min(numel(ydata),posy(end)+offset);

	data  = double(imread(remapath,'PixelRegion',{[id1y,id2y],[id1x,id2x]}));
	xdata=xdata(id1x:id2x);
	ydata=ydata(id1y:id2y);
end

%convert no coverage data
data(find(data==-9999))=NaN;

sout = InterpFromGrid(xdata,ydata,data,X,Y);
