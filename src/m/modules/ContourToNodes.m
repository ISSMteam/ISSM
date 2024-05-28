function flags = ContourToNodes(x,y,contourname,edgevalue);
%CONTOURTONODES - flags vertices inside contour
%
%   Usage:
%      flags = ContourToNodes(x,y,contourname,edgevalue);
%
%   x,y: list of nodes
%   contourname: name of Argus or Shape file containing the contours, or resulting structure from call to expread
%   edgevalue: integer (0, 1 or 2) defining the value associated to the nodes on the edges of the polygons
%   flags: vector of flags (0 or 1), of size nodes

%Check usage
if nargin~=4,
	help ContourToNodes
	error('Wrong usage (see above)');
end

%Some conversion of files:  
if ischar(contourname),
	[path,name,ext]=fileparts(contourname); 
	if strcmpi(ext,'.shp'),
		%read contour from shapefile:
		contourname=shpread(contourname); 

		%FIXME: I don't think we need to convert it to a file, ContourToMesh_matlab should be able to take a structure
		%write it to a temporary filename: 
		%contourname=[tempname '.exp'];
		%expwrite(contour,contourname);
	end
end

%Call mex module
[flags] = ContourToNodes_matlab(x,y,contourname,edgevalue);
