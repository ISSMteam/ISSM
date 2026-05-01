function [dfdx, dfdy, df] = slope(md, f)
%SLOPE - compute the gradient of any field
%
%   Usage:
%      [dfdx, dfdy, df] = slope(md)
%      [dfdx, dfdy, df] = slope(md, md.results.TransientSolution(1).Surface)
%                    df = slope(md, md.geometry.surface);
%
%  where df = sqrt(dfdx^2 + dfdy^2)

%load some variables (it is much faster if the variab;es are loaded from md once for all) 
if dimension(md.mesh)==2
	numberofelements = md.mesh.numberofelements;
	numberofnodes    = md.mesh.numberofvertices;
	index            = md.mesh.elements;
	x                = md.mesh.x;
   y                = md.mesh.y;
else
	numberofelements = md.mesh.numberofelements2d;
	numberofnodes    = md.mesh.numberofvertices2d;
	index            = md.mesh.elements2d;
	x                = md.mesh.x2d;
	y                = md.mesh.y2d;
end

if nargin==1
	f = md.geometry.surface;
end

%compute nodal functions coefficients N(x,y)=alpha x + beta y + gamma
[alpha beta]=GetNodalFunctionsCoeff(index,x,y);

summation=[1;1;1];
dfdx = (f(index).*alpha)*summation;
dfdy = (f(index).*beta)*summation;
if dimension(md.mesh)==3
	dfdx = project3d(md,'vector',dfdx,'type','element');
	dfdy = project3d(md,'vector',dfdy,'type','element');
end

%Compute magnitude
df   = sqrt(dfdx.^2+dfdy.^2);

%return magnitude only
if nargout==1; dfdx = df; end
