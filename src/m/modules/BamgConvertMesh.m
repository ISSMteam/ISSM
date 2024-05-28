function [bamggeom, bamgmesh] = BamgConvertMesh(index,x,y);
%BAMGCONVERTMESH - Convert [x y index] to a bamg geom and mesh geom
%   
%   Usage:
%      [bamggeom, bamgmesh] = BamgConvertMesh(index,x,y);
%   
%   index: index of the mesh
%   x,y: coordinates of the nodes

% Check usage
if nargin~=3
	help BamgConvertMesh
	error('Wrong usage (see above)');
end

% Call mex module
[bamggeom, bamgmesh] = BamgConvertMesh_matlab(index,x,y);
