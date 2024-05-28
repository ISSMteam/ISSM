function [output] = interpRTopo2(X,Y,varargin),
%INTERPRTOPO2 - interp from RTOPO-2 onto X and Y
%
%   Usage:
%      bed = interpRTopo2(X,Y,varargin),
%
%   varargin = 1 (Greenland), default
%             -1 (Antarctica)

switch oshostname(),
	case {'ronne'}
		rootname='/home/ModelData/Global/RTopo-2/RTopo-2.0.1_30sec_bedrock_topography.nc';
	case {'totten'}
		rootname='/totten_1/ModelData/Global/RTopo-2/RTopo-2.0.1_30sec_bedrock_topography.nc';
	otherwise
		error('machine not supported yet');
end
verbose = 1;

if nargin==3,
	hemisphere = varargin{1};
else
	hemisphere = +1;
end
if abs(hemisphere)~=1,
	error('hemisphere should be +/-1');
end

if hemisphere==+1,
	if verbose, disp('   -- RTopo-2: convert to lat/lon using Greenland projection'); end
	[LAT,  LON  ] = xy2ll(double(X(:)),double(Y(:)),+1,45,70);
else
	if verbose, disp('   -- RTopo-2: convert to lat/lon using Antarctica projection'); end
	[LAT,  LON  ] = xy2ll(double(X(:)),double(Y(:)),-1,0,71);
end

Y=reshape(LAT,size(X)); X=reshape(LON,size(X));

xdata = double(ncread(rootname,'lon'));
ydata = double(ncread(rootname,'lat'));

offset=2;

xmin=min(X(:)); xmax=max(X(:));
posx=find(xdata<=xmax);
id1x=max(1,find(xdata>=xmin,1)-offset);
id2x=min(numel(xdata),posx(end)+offset);

ymin=min(Y(:)); ymax=max(Y(:));
posy=find(ydata<=ymax);
id1y=max(1,find(ydata>=ymin,1)-offset);
id2y=min(numel(ydata),posy(end)+offset);

if verbose, disp('   -- RTopo-2: reading bed topography'); end
data  = double(ncread(rootname,'bedrock_topography',[id1x id1y],[id2x-id1x+1 id2y-id1y+1],[1 1]))';
xdata=xdata(id1x:id2x);
ydata=ydata(id1y:id2y);
data(find(data==-9999))=NaN;

if verbose, disp('   -- RTopo-2: interpolating'); end
output = InterpFromGrid(xdata,ydata,data,double(X),double(Y));
