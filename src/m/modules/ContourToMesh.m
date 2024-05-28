function varargout = ContourToMesh(index,x,y,contourname,interptype,edgevalue);
%CONTOURTOMESH - Flag the elements or nodes inside a contour
%
%   Usage:
%      [in_nod,in_elem]=ContourToMesh(index,x,y,contourname,interptype,edgevalue);
%	
%   index,x,y: mesh triangulation
%   contourname: name of .exp or .shp file containing the contours.
%   interptype: string defining type of interpolation ('element', or 'node').
%   edgevalue: integer (0, 1, or 2) defining the value associated to the nodes on the edges of the polygons.
%   in_nod: vector of flags (0 or 1), of size nods if interptype is set to 'node' or 'element and node',
%      or of size 0 otherwise
%   in_elem: vector of flags (0 or 1), of size nel if interptype is set to 'element' or 'element and node',
%      or of size 0 otherwise.    
%
%   Example:
%      in_nod=ContourToMesh(md.elements,md.x,md.y,'Contour.exp','node',1)
%      in_elements=ContourToMesh(md.elements,md.x,md.y,'Contour.exp','element',0)
%      [in_nodes,in_elements]=ContourToMesh(md.elements,md.x,md.y,'Contour.exp','element and node',0)

%Check usage
if nargin~=6
	help ContourToMesh
	error('Wrong usage(see above)');
end

%Some conversion of files: 
if ischar(contourname),
	[path,name,ext]=fileparts(contourname); 
	if strcmpi(ext,'.shp'),
		%read contour from shapefile
		contourname=shpread(contourname); 

		%FIXME: I don't think we need to convert it to a file, ContourToMesh_matlab should be able to take a structure
		%write it to a temporary filename: 
		%contourname=[tempname '.exp'];
		%expwrite(contour,contourname);
	end
end

%Call mex module
[in_nod,in_elem] = ContourToMesh_matlab(index,x,y,contourname,interptype,edgevalue);

switch(interptype)
	case 'element'
		varargout{1} = in_elem;
	case 'node'
		varargout{1} = in_nod;
	case 'node and element'
		varargout{1} = in_nod;
		varargout{2} = in_elem;
	otherwise
		error(['interpolation type ''' interptype ''' not supported yet']);
end
