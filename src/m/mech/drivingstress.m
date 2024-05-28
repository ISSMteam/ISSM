function [px,py,pmag]=drivingstress(md)
%DRIVINGSTRESS -  evaluates the driving stress
%
%   The driving stress is computed according to the following formula: 
%   driving stress= rho_ice*g*H*slope
%
%   Usage:
%      [Fx,Fy,Fmag]=drivingstress(md)

%Get slope
[sx,sy,s]=slope(md);

%Average thickness over elements
thickness_bar=(md.geometry.thickness(md.mesh.elements(:,1))+md.geometry.thickness(md.mesh.elements(:,2))+md.geometry.thickness(md.mesh.elements(:,3)))/3;

px=-md.materials.rho_ice*md.constants.g*thickness_bar.*sx;
py=-md.materials.rho_ice*md.constants.g*thickness_bar.*sy;
pmag=sqrt(px.^2+py.^2);
