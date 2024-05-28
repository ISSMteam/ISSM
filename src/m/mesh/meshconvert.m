function md=meshconvert(md,varargin)
%CONVERTMESH - convert mesh to bamg mesh
%
%   Usage:
%      md=meshconvert(md);
%      md=meshconvert(md,index,x,y);

if nargin~=1 & nargin~=4,
	help meshconvert
	error('meshconvert error message: bad usage');
end

if nargin==1,
	index = md.mesh.elements;
	x     = md.mesh.x;
	y     = md.mesh.y;
else
	index = varargin{1};
	x     = varargin{2};
	y     = varargin{3};
end

%call Bamg
[bamgmesh_out bamggeom_out]=BamgConvertMesh(index,x,y);

% plug results onto model
md.private.bamg          = struct();
md.private.bamg.mesh     = bamgmesh(bamgmesh_out);
md.private.bamg.geometry = bamggeom(bamggeom_out);
md.mesh.x              = bamgmesh_out.Vertices(:,1);
md.mesh.y              = bamgmesh_out.Vertices(:,2);
md.mesh.elements       = bamgmesh_out.Triangles(:,1:3);
md.mesh.edges          = bamgmesh_out.IssmEdges;
md.mesh.segments       = bamgmesh_out.IssmSegments(:,1:3);
md.mesh.segmentmarkers = bamgmesh_out.IssmSegments(:,4);

%Fill in rest of fields:
md.mesh.numberofelements = size(md.mesh.elements,1);
md.mesh.numberofvertices = length(md.mesh.x);
md.mesh.numberofedges    = size(md.mesh.edges,1);
md.mesh.vertexonboundary = zeros(md.mesh.numberofvertices,1); md.mesh.vertexonboundary(md.mesh.segments(:,1:2)) = 1;
md.mesh.elementconnectivity=md.private.bamg.mesh.ElementConnectivity;
md.mesh.elementconnectivity(find(isnan(md.mesh.elementconnectivity)))=0;
