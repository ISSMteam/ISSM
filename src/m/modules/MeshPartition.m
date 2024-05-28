function [element_partitioning, node_partitioning] = MeshPartition(md,numpartitions)
%MESHPARTITION - Partition mesh according to the number of areas, using Metis library.
%
%	   Usage:
%			[element_partitioning,node_partitioning]=MeshPartition(md,numpartitions);
%
%	   element_partitioning: Vector of partitioning area numbers, for every element.
%	   node_partitioning: Vector of partitioning area numbers, for every node.

% Check usage
if nargin~=2
	help MeshPartition
	error('Wrong usage (see above)');
end

%Get mesh info from md.mesh
numberofvertices = md.mesh.numberofvertices;
numberofelements = md.mesh.numberofelements;
elements         = md.mesh.elements;
numberofelements2d = 0;
numberofvertices2d = 0;
numberoflayers     = 1;
elements2d         = [];
if isa(md.mesh,'mesh3dprisms')
	elementtype = 'Penta';
	numberofelements2d = md.mesh.numberofelements2d;
	numberofvertices2d = md.mesh.numberofvertices2d;
	numberoflayers     = md.mesh.numberoflayers;
	elements2d         = md.mesh.elements2d;
elseif isa(md.mesh,'mesh2d')
	elementtype = 'Tria';
elseif isa(md.mesh,'mesh2dvertical')
	elementtype = 'Tria';
end

% Call mex module
[element_partitioning, node_partitioning] = MeshPartition_matlab(numberofvertices,elements,numberofvertices2d,elements2d,numberoflayers,elementtype,numpartitions);
