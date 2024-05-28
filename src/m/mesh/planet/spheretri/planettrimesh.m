function md=planettrimesh(md,shape,radius,refinement)
%PLANETTRIMESH: build 2d shell mesh
%
% Usage: md=planettrimesh(md,shape,radius,refinement)
%

results = sphere_tri(shape,refinement,radius);
md.mesh=mesh3dsurface(); %???
md.mesh.x=results.vertices(:,1);
md.mesh.y=results.vertices(:,2);
md.mesh.z=results.vertices(:,3);
md.mesh.elements=results.faces;

md.mesh.r=sqrt(md.mesh.x.^2+md.mesh.y.^2+md.mesh.z.^2);
md.mesh.lat=acos(md.mesh.z./md.mesh.r);
md.mesh.long=atan2(md.mesh.y,md.mesh.x);

md.mesh.numberofvertices=length(md.mesh.x);
md.mesh.numberofelements=size(md.mesh.elements,1);
