function showbasins(varargin)
%SHOWBASINS - return basins that are within the xlim and ylim
%
%   Usage:
%      names=showbasins(options);
%   Options: 
%      'unit' default 1
%      'hemisphere': default +1;
%      'central_meridian: 45 for Greenland and 0 for Antarctica
%      'standard_parallel: 70 for Greenland and 71 for Antarctica
%

%is varargin an options database already?
if nargin==0,
	options=pairoptions(varargin{:});
elseif (isa(varargin{1},'plotoptions') | isa(varargin{1},'pairoptions')),
	%do nothing to the varargin: 
	options=varargin{1};
else
	%process varargin for options: 
	options=pairoptions(varargin{:});
end

%recover some options, and set defaults
unitmultiplier=getfieldvalue(options,'unit',1);
fontsize=getfieldvalue(options,'fontsize',12);
hemisphere=getfieldvalue(options,'hemisphere');

if strcmpi(hemisphere,'s'),
	hemisphere=-1;
elseif strcmpi(hemisphere,'n'),
	hemisphere=+1;
else
	error('showbasins error message: hemispehre should be either ''n'' or ''s''');
	end

if hemisphere==+1,
	central_meridian=getfieldvalue(options,'central_meridian',45);
	standard_parallel=getfieldvalue(options,'standard_parallel',70);
else
	central_meridian=getfieldvalue(options,'central_meridian',0);
	standard_parallel=getfieldvalue(options,'standard_parallel',71);
end

%Ok, find basin we are talking about: 
load([jplsvn '/projects/ModelData/Names/Names.mat']);

%Get xlim and ylim, and convert into lat,long: 
xlimits=xlim; x0=xlimits(1); x1=xlimits(2);
ylimits=ylim; y0=ylimits(1); y1=ylimits(2);

%Convert names lat and long into x,y:
lat=cell2mat(names(:,3));
long=cell2mat(names(:,2));

%Now, convert lat,long into x,y:
[x,y]=ll2xy(lat,long,hemisphere,central_meridian,standard_parallel);

%Find  x,y within xlimits and ylimits: 
locations=find(x>x0 & x<x1 & y>y0 & y<y1);

%Go through locations, and display the names: 
for i=1:size(locations,1),
	hold on,
	plot(x(locations(i)),y(locations(i)),'r.');
	t=text(x(locations(i)),y(locations(i)),names{locations(i),1}); 
	set(t,'FontSize',fontsize);
end
