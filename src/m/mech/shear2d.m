function [sx,sy,sxy,s]=shear2d(md)
%SHEAR2D - computes 2d strain rate
%
%   This routine computes the strain rate of 2d models
%
%   Usage:
%      [sx,sy,sxy,s]=shear2d(md);
%      s=shear2d(md);

[alpha beta]=GetNodalFunctionsCoeff(md.mesh.elements,md.mesh.x,md.mesh.y); 

summation=[1;1;1];
sx=(md.initialization.vx(md.mesh.elements).*alpha)*summation;
uy=(md.initialization.vx(md.mesh.elements).*beta)*summation;
vx=(md.initialization.vy(md.mesh.elements).*alpha)*summation;
sy=(md.initialization.vy(md.mesh.elements).*beta)*summation;						
sxy=(uy+vx)/2;
s=sqrt(sx.^2+sy.^2+sxy.^2+sx.*sy);

%if user requested only one output, it must be the norm
if nargout==1,
	sx=s;
end
