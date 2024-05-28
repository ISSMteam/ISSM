function md=adjacency(md)
%ADJACENCY -  compute adjacency matrix, list of vertices and list of weights.
%
%  function to create the adjacency matrix from the connectivity table.
%
%  the required output is:
%    md.adj_mat     (double [sparse nv x nv], vertex adjacency matrix)
%    md.qmu.vertex_weight        (double [nv], vertex weights)

indi=[md.mesh.elements(:,1);md.mesh.elements(:,2);md.mesh.elements(:,3)];
indj=[md.mesh.elements(:,2);md.mesh.elements(:,3);md.mesh.elements(:,1)];
values=1;

md.qmu.adjacency=sparse(indi,indj,values,md.mesh.numberofvertices,md.mesh.numberofvertices);
md.qmu.adjacency=double([md.qmu.adjacency | md.qmu.adjacency']);

%now, build vwgt:
areas=GetAreas(md.mesh.elements,md.mesh.x,md.mesh.y);

%get node connectivity
md.mesh.vertexconnectivity=NodeConnectivity(md.mesh.elements,md.mesh.numberofvertices);

connectivity=md.mesh.vertexconnectivity(:,1:end-1);
pos=find(connectivity);
connectivity(pos)=areas(connectivity(pos))/3;
md.qmu.vertex_weight=sum(connectivity,2);
