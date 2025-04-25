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
if isfield(Tinfo(1),'ModelPixelScaleTag')
	dx    = Tinfo(1).ModelPixelScaleTag(1);
	dy    = Tinfo(1).ModelPixelScaleTag(2);
	minx  = Tinfo(1).ModelTiepointTag(4);
	maxy  = Tinfo(1).ModelTiepointTag(5);
	assert(dx>0); assert(dy>0);
elseif isfield(Tinfo(1),'ModelTransformationTag')
	dx   = Tinfo(1).ModelTransformationTag(1);
	dy   = -Tinfo(1).ModelTransformationTag(6);
	minx = Tinfo(1).ModelTransformationTag(4);
	maxy = Tinfo(1).ModelTransformationTag(8);
	assert(dx>0); assert(dy<0);
else
   error('image info cannot be retrieved for this geotiff');
end

%Generate vectors
xm = minx + dx/2 + ((0:N-1).*dx);
ym = maxy - dy/2 - ((M  -1:-1:0).*dy);

%Read image
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
