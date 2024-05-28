function labels = labelconnectedregions(md)
%LABELCONNECTEDREGIONS - label connected components of a mesh
%
%   Usage:
%      labels = labelconnectedregions(md)

if size(md.mesh.elements,2)~=3,
	error('not suppored yet (but easy to extend :)');
end

disp('Generate adjacency matrix');
pairs = [
md.mesh.elements(:,[1 2])
md.mesh.elements(:,[2 1])
md.mesh.elements(:,[2 3])
md.mesh.elements(:,[3 2])
md.mesh.elements(:,[3 1])
md.mesh.elements(:,[1 3])
];
A = sparse(pairs(:,1), pairs(:,2), 1);

disp('Construct graph');
G = graph(A);

disp('Label connected pieces');
labels = conncomp(G);
