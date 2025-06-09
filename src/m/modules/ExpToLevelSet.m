function distance = ExpToLevelSet(x,y,contourname,dateformat);
%EXPTOLEVELSET - Determine levelset distance between a contour and a cloud of points
%
%   Usage:
%      distance=ExpToLevelSet(x,y,contourname);
%
%   x,y:	cloud point
%   contourname:	name of .exp file containing the contours
%	 dateformat: not in use, by default
%   distance:	distance vector representing a levelset where the 0 level is one of the contour segments
%
%   Example:
%      distance=ExpToLevelSet(md.mesh.x,md.mesh.y,'Contour.exp');

% Check usage
if nargin~=3 & nargin~=4
	help ExpToLevelSet
	error('Wrong usage (see above)');
end
if nargin < 4
	dateformat = '';
end

multipleShp = 0;

if ischar(contourname),
	[path,name,ext]=fileparts(contourname);
	if strcmpi(ext,'.shp'),
		%read contour from shapefile
		contourname=shpread(contourname);
		if isstruct(contourname)
			multipleShp = size(contourname, 2);
		end
	end
end

% Call mex module
if multipleShp>0 
	% shp file contains multiple contours
	distance = zeros(length(x)+1, multipleShp);
	for i = 1:multipleShp
		distance(1:end-1, i) = ExpToLevelSet_matlab(x,y,contourname(i));
		% append the NAME information at the end of the distance
		if isfield(contourname(i),'NAME')
			distance(end, i) = contourname(i).NAME;
		elseif isfield(contourname(i),'Date')
			if strcmp(dateformat, '')
				distance(end, i) = date2decyear(datenum(contourname(i).Date));
			else
				distance(end, i) = date2decyear(datenum(contourname(i).Date, dateformat));
			end
			% to deal with 
		else
			distance(end, i) = i;
		end        
	end
else
	% single shape or exp file
	distance = ExpToLevelSet_matlab(x,y,contourname);
end
