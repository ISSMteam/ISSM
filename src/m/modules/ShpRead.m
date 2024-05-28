function structout = ShpRead(filename);
%	SHPREAD - Read shapefile
%	
%	   This module reads shapefiles and converts them to matlab/python structures
%	
%	   Usage:
%	      structout = ShpRead(filename);
%	
%	   Examples:
%	      structout = ShpRead('file.shp');

% Check usage
if nargin~=1
	help ShpRead
	error('Wrong usage: No file specified');
end

% Call mex module
structout = ShpRead_matlab(filename);
