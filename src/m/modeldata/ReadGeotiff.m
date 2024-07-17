function [xm ym data] = ReadGeotiff(geotiffname,nanValue,fillholes)
%READGEOTIFF - read geotiff and output coordinates and map
%
%   Usage:
%      [xm ym data] = ReadGeotiff(geotiffname);
%      [xm ym data] = ReadGeotiff(geotiffname,nanValue,fillholes)


if nargin < 2
	nanValue = 10^30;
	fillholes = false;
end
if nargin < 3
	fillholes = false;
end

%Get image info
Tinfo = imfinfo(geotiffname);
N     = Tinfo(1).Width;
M     = Tinfo(1).Height;
dx    = Tinfo(1).ModelPixelScaleTag(1);
dy    = Tinfo(1).ModelPixelScaleTag(2);
minx  = Tinfo(1).ModelTiepointTag(4);
maxy  = Tinfo(1).ModelTiepointTag(5);

%Generate vectors
xm = minx + dx/2 + ((0:N-1).*dx);
ym = maxy - dy/2 - ((M  -1:-1:0).*dy);

%Read image
assert(dx>0); assert(dy>0);
ym = fliplr(ym);

data  = double(imread(geotiffname));
if nanValue > 0
	data(find(abs(data)>=nanValue))=NaN;
else 
	data(find(data<=nanValue))=NaN;
end
if fillholes
	disp('Filling holes');
	data = inpaint_nans(data);
	disp('done');
end

data(data==-9999)=NaN;
