function varargout=basinzoom(varargin)
%ANTZOOM - zoom on a basin in Antarctica or Greenland.
%
%   This function zooms on an existing figure describing Antarctica or Greenland 
%   The zooming depends on the region name provided as input. 
%
%   Usage:
%      varargout=basinzoom(options)

%recover some options, and set defaults

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

unitmultiplier=getfieldvalue(options,'unit',NaN);
basin=getfieldvalue(options,'basin');

if exist(options,'basindelta'),

	basindeltax=getfieldvalue(options,'basindelta',300); 
	basindeltay=getfieldvalue(options,'basindelta',300); 
else
	basindeltax=getfieldvalue(options,'basindeltax',300); 
	basindeltay=getfieldvalue(options,'basindeltay',300);
end

%multiply by 1000 to get kms
basindeltax=basindeltax*1000;
basindeltay=basindeltay*1000;

%Ok, find basin we are talking about: 
load([jplsvn() '/ModelData/Names/Names.mat']);

%Go through names: 
found=0;
for i=1:size(names,1),
	if strcmpi(names{i,1},basin),
		%ok, we've got the region. Get lat and long: 
		long=names{i,2};
		lat=names{i,3};
		hemisphere=names{i,4};
		found=1;
		break;
	end
end

if ~found,
	error(['basinzoom error message: cannot find basin ' basin '. Use isbasin to determine a basin name.']);
end

if hemisphere==+1,
	central_meridian=getfieldvalue(options,'central_meridian',45);
	standard_parallel=getfieldvalue(options,'standard_parallel',70);
else
	central_meridian=getfieldvalue(options,'central_meridian',0);
	standard_parallel=getfieldvalue(options,'standard_parallel',71);
end

%Transform lat long into x,y: 
[xc,yc]=ll2xy(lat,long,hemisphere,central_meridian,standard_parallel);

%compute x0,x1 and y0,y1 using basindeltax and basindeltay
x0=xc-basindeltax/2;
x1=xc+basindeltax/2;
y0=yc-basindeltay/2;
y1=yc+basindeltay/2;

if ~isnan(unitmultiplier)
	x0=x0*unitmultiplier;
	x1=x1*unitmultiplier;
	y0=y0*unitmultiplier;
	y1=y1*unitmultiplier;
end

%if output arguments are present, return the limits, 
%otherwise, set them on the current graphic. 
if nargout==2,
	found=1;
	varargout{1}=[x0 x1];
	varargout{2}=[y0 y1];
else
	xlim([x0 x1]);
	ylim([y0 y1]);
	found=1;
	daspect([1;1;1]);
end
