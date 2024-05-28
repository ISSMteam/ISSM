function div=divergence(md,a,b)
%DIVERGENCE - divergence of [a;b] vector, using model's triangulation.
%
%   Usage:
%      div=divergence(md,a,b)

if (dimension(md.mesh)==2),
	numberofelements=md.mesh.numberofelements;
	numberofnodes=md.mesh.numberofvertices;
	index=md.mesh.elements;
	x=md.mesh.x; y=md.mesh.y; z=md.mesh.z;
else
	numberofelements=md.mesh.numberofelements2d;
	numberofnodes=md.mesh.numberofvertices2d;
	index=md.mesh.elements2d;
	x=md.mesh.x2d; y=md.mesh.y2d;
end

%compute nodal functions coefficients N(x,y)=alpha x + beta y + gamma
[alpha beta]=GetNodalFunctionsCoeff(index,x,y);

summation=[1;1;1];
dx=(a(index).*alpha)*summation;
dy=(b(index).*beta)*summation;

div=dx+dy;
div=averaging(md,div,1);
