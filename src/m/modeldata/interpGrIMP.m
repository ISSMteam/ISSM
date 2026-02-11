function sout = interpGrIMP(X,Y)

switch oshostname(),
	case {'totten'}
		datapath='/totten_1/ModelData/Greenland/GrIMP/GrIMP_100m.tif'; %McGregor
		datapath='/totten_1/ModelData/Greenland/GrIMP/GIMP2_3413.tif'; %Bea
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
	[data,R] = geotiffread(datapath);
	data=double(flipud(data));
	xdata=R.XLimWorld(1):R.DeltaX:R.XLimWorld(2); xdata=xdata(:);
	xdata =(xdata(1:end-1)+xdata(2:end))/2;
	ydata=R.YLimWorld(2):R.DeltaY:R.YLimWorld(1); ydata=flipud(ydata(:));
	ydata =(ydata(1:end-1)+ydata(2:end))/2;
else

	%Get image info
	Tinfo = imfinfo(datapath);
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

	data  = double(imread(datapath,'PixelRegion',{[id1y,id2y],[id1x,id2x]}));
	xdata=xdata(id1x:id2x);
	ydata=ydata(id1y:id2y);
end

%corrections
data(data==-9999) = NaN;
data(data>1e6) = NaN;
sout = InterpFromGrid(xdata,ydata,data,X,Y,'linear');
sout(sout==-9999) = NaN;
