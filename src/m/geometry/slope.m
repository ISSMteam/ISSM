function [sx,sy,s]=slope(md,surf)
%SLOPE - compute the surface slope
%
%   Usage:
%      [sx,sy,s]=slope(md)
%      [sx,sy,s]=slope(md,md.results.TransientSolution(1).Surface)

%load some variables (it is much faster if the variab;es are loaded from md once for all) 
if dimension(md.mesh)==2,
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

if nargin==1,
	surf=md.geometry.surface;
end
%compute nodal functions coefficients N(x,y)=alpha x + beta y + gamma
[alpha beta]=GetNodalFunctionsCoeff(index,x,y);

summation=[1;1;1];
sx=(surf(index).*alpha)*summation;
sy=(surf(index).*beta)*summation;
s=sqrt(sx.^2+sy.^2);

if dimension(md.mesh)==3,
	sx=project3d(md,'vector',sx,'type','element');
	sy=project3d(md,'vector',sy,'type','element');
	s=sqrt(sx.^2+sy.^2);
end
