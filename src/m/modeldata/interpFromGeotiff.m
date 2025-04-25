function dataout = interpFromGeotiff(geotiffname, X, Y, nanValue, fillholes, method)
%INTERPFROMGEOTIFF - interpolate field in geotiff onto list of points
%
%   Usage:
%      dataout = interpFromGeotiff(geotiffname,X,Y,nanValue,fillholes,method)
%      dataout = interpFromGeotiff(geotiffname,X,Y);
%
%   Options:
%      nanValue:  which values are considered NaN? default is 1e30
%      fillholes: default is false
%      method:    default is 'linear'

%Default options
if nargin < 6
	method = 'linear';
end
if nargin < 5
	fillholes = false;
end
if nargin < 4
	nanValue = 10^30;
end


usemap = 0;
if license('test','map_toolbox')==0
	disp('WARNING: map toolbox not installed, trying in-house code');
	usemap = 0;
elseif license('checkout','map_toolbox')==0
	disp('WARNING: map toolbox not available (checkout failed), trying in-house code');
	usemap = 0;
end

if usemap
	[data,R] = geotiffread(geotiffname);
	data=double(flipud(data));
	xdata=R.XLimWorld(1):R.DeltaX:R.XLimWorld(2); xdata=xdata(:);
	xdata =(xdata(1:end-1)+xdata(2:end))/2;
	ydata=R.YLimWorld(2):R.DeltaY:R.YLimWorld(1); ydata=flipud(ydata(:));
	ydata =(ydata(1:end-1)+ydata(2:end))/2;

else

	%Get image info
	Tinfo = imfinfo(geotiffname);
	N     = Tinfo(1).Width;
	M     = Tinfo(1).Height;
	dx    = Tinfo(1).ModelPixelScaleTag(1);
	dy    = Tinfo(1).ModelPixelScaleTag(2);
	minx  = Tinfo(1).ModelTiepointTag(4);
	maxy  = Tinfo(1).ModelTiepointTag(5);

	%Generate vectors
	xdata = minx + dx/2 + ((0:N-1).*dx);
	ydata = maxy - dy/2 - ((M  -1:-1:0).*dy);

	%Read image
	if 1
		assert(dx>0); assert(dy>0);
		ydata = fliplr(ydata);

		%Get pixels we are interested in
		offset=2;

		xmin=min(X(:)); xmax=max(X(:));
		xflags = (xdata>xmin & xdata<xmax);
		if ~any(xflags); dataout = NaN(size(X)); return; end
		id1x=find(xflags, 1, 'first'); id1x = max(1,id1x -offset);
		id2x=find(xflags, 1, 'last' ); id2x = min(numel(xdata),id2x+offset);

		ymin=min(Y(:)); ymax=max(Y(:));
		yflags = (ydata>ymin & ydata<ymax);
		if ~any(yflags); dataout = NaN(size(X)); return; end
		id1y=find(yflags, 1, 'first'); id1y = max(1,id1y -offset);
		id2y=find(yflags, 1, 'last' ); id2y = min(numel(ydata),id2y+offset);

		data  = double(imread(geotiffname,'PixelRegion',{[id1y,id2y],[id1x,id2x]}));
		xdata = xdata(id1x:id2x);
		ydata = ydata(id1y:id2y);
	else
		data=double(flipud(imread(geotiffname)));
	end
end

%Deal with NaNs
if nanValue > 0
	data(find(abs(data)>=nanValue))=NaN;
else 
	data(find(data<=nanValue))=NaN;
end

%Fill holes using inpaint_nans
if fillholes
	disp('Filling holes');
	data = inpaint_nans(data);
	disp('done');
end

%Interpolate
dataout = InterpFromGrid(xdata, ydata, data, X, Y, method);
dataout(dataout==-9999)=NaN;
