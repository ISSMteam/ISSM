function [vx, vy, vz, vel]=sia(md)
%SIA - computation of Shallow Ice velocities
%
%   This routine uses the model of SIA to compute the velocities
%   of a 2d model using the surface slope
%
%   Usage:
%      [vx, vy, vz, vel]=sia(md)

disp('Info: assuming no basal sliding and isothermal ice');

n=3;

rhog = (md.materials.rho_ice*md.constants.g);

if md.mesh.dimension==2,
	%Get slope
	[sx,sy,s]=slope(md);

	%Average thickness and B over all elements.
	summer=[1;1;1];
	hel=md.geometry.thickness(md.mesh.elements)*summer/3;
	Bel=md.materials.rheology_B(md.mesh.elements)*summer/3;

	Ael=Bel.^(-3);

	vx=-2*rhog^n*s.^(n-1).*sx.*Ael/(n+2).*hel.^(n+1);
	vy=-2*rhog^n*s.^(n-1).*sy.*Ael/(n+2).*hel.^(n+1);
	vz = zeros(size(vx));
else

	%Smooth surface slope a little bit first
	[sx,sy,s]=slope(md);
	sx = averaging(md, sx, 3);
	sy = averaging(md, sy, 3);
	surf = sqrt(sx.^2 + sy.^2);

	%Get other variables
	s  = md.geometry.surface;
	z  = md.mesh.z;
	H  = md.geometry.thickness;

	A = md.materials.rheology_B.^-3;

	vx =-2*rhog^n*surf.^(n-1).*sx.*A/(n+1).*(H.^(n+1) - (s-z).^(n+1));
	vy =-2*rhog^n*surf.^(n-1).*sy.*A/(n+1).*(H.^(n+1) - (s-z).^(n+1));

	chi = 1-(s-z)./H;
	vz  = -md.smb.mass_balance.*((n+2)*chi + (1-chi).^(n+2) -1)/(n+1);
end

vx = vx*md.constants.yts;
vy = vy*md.constants.yts;

vel=sqrt(vx.^2+vy.^2+vz.^2);
