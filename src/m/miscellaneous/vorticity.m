function rot = vorticity(md,vx,vy)
%VORTICITY - calculates 2d vorticity
%
%   rot = d/dx(vy) - d/dy(vx)
%
%   Usage:
%      rot = vorticity(vx,vy)


%load some variables (it is much faster if the variab;es are loaded from md once for all) 
if ~strcmpi(md.mesh.domaintype(),'3D'),
	numberofelements=md.mesh.numberofelements;
	numberofnodes=md.mesh.numberofvertices;
	index=md.mesh.elements;
	x=md.mesh.x; y=md.mesh.y;
else
	numberofelements=md.mesh.numberofelements2d;
	numberofnodes=md.mesh.numberofvertices2d;
	index=md.mesh.elements2d;
	x=md.mesh.x2d; y=md.mesh.y2d;
end

%compute nodal functions coefficients N(x,y)=alpha x + beta y + gamma
[alpha beta]=GetNodalFunctionsCoeff(index,x,y);

summation=[1;1;1];
dvydx=(vy(index).*alpha)*summation;
dvxdy=(vx(index).*beta)*summation;
rot=dvxdy - dvydx;

if strcmpi(domaintype(md.mesh),'3D'),
	rot=project3d(md,'vector',rot,'type','element');
end
