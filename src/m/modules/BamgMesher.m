function [bamgmesh,bamggeom] = BamgMesher(bamgmesh,bamggeom,bamgoptions);
%BAMGMESHER
%
%   Usage:
%      [bamgmesh, bamggeom] = BamgMesher(bamgmesh,bamggeom,bamgoptions);
%
%   bamgmesh: input bamg mesh
%   bamggeom: input bamg geometry for the mesh
%   bamgoptions: options for the bamg mesh

% Check usage
if nargin~=3
	help BamgMesher
	error('Wrong usage (see above)');
end

% Call mex module
[bamgmesh, bamggeom] = BamgMesher_matlab(bamgmesh,bamggeom,bamgoptions);

