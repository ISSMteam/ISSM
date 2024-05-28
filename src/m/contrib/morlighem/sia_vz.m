function [velz]=sia_vz(md)
%SIA_VZ - computation vertical speed based on Shallow Ice Approximation
%
%   Usage:
%      [vz]=sia_vz(md)

if md.mesh.dimension~=3
	error('Only 3d meshes are allowed to compute vz');
end

ws = - md.surfaceforcings.mass_balance;
n  = md.materials.rheology_n(1); %just take the first one
z  = md.mesh.z;
b  = md.geometry.base;
H  = md.geometry.thickness;
s  = md.geometry.surface;

vz = ws.*(n+2)/(n+1).*((z-b)./H + 1/(n+2)*(((s-z)./H).^(n+2)-1));
