function [segments] = MeshProfileIntersection(index,x,y,filename);
%MESHPROFILEINTERSECTION - Take a .exp file (made of several profiles), and figures out its intersection with a mesh.
%
%	   usage:
%	   [segments]=MeshProfileIntersection(index,x,y,filename);
%
%	   input:
%	        index,x,y is a triangulation
%	        filename: name of Argus or Shape file containing the segments (can be groups of disconnected segments)
%	   output:
%	        segments: array made of x1,y1,x2,y2,element_id lines (x1,y1) and (x2,y2) are segment extremities for a segment 
%	        belonging to the elemnt_id element. there are as many lines in segments as there are segments intersecting the 
%	        mesh.

% Check usage
if nargin~=4
	help MeshProfileIntersection
	error('Wrong usage (see above)');
end

[path,name,ext]=fileparts(filename); 
if strcmpi(ext,'.shp'),
	
	%convert to expfile and store in a temporary directory: 
	oldfilename=filename;
	filename=[tempname '.exp'];
	shp2exp(oldfilename,filename);
end

% Call mex module
[segments] = MeshProfileIntersection_matlab(index,x,y,filename);
