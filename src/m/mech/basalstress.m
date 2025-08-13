function [bx by b]=basalstress(md)
%BASALSTRESS - compute basal stress from basal drag and geometric information. 
%
%      Computes basal stress from geometric information and ice velocity in md.initialization. Follows the basal stress definition in "src/c/classes/Loads/Friction.cpp", lines 1102-1136.
%
%   Usage:
%      [bx by b]=basalstress(md);
%
%   See also: plot_basaldrag

%Check md.friction class
if ~isa(md.friction,'friction')
	error('Error: md.friction only supports "friction.m" class.');
end

%compute exponents
s=averaging(md,1./md.friction.p,0);
r=averaging(md,md.friction.q./md.friction.p,0);

%Compute effective pressure
g        =md.constants.g;
rho_ice  =md.materials.rho_ice;
rho_water=md.materials.rho_water;

sealevel=0;
p_ice=g*rho_ice*md.geometry.thickness;
switch(md.friction.coupling)
   case 0
		p_water=g*rho_water*(sealevel-md.geometry.base);
      N = p_ice-p_water;
   case 1
      N = p_ice;
   case 2
		p_water=g*rho_water*(sealevel-md.geometry.base);
		p_water=max(p_water,0.0);
		N = p_ice-p_water;
   case 3
      N = max(md.friction.effective_pressure, 0);
   case 4
      error('md.friction.coupling=4 is not supported yet.');
   otherwise
      error('not supported yet');
end

%compute sliding velocity
ub=sqrt(md.initialization.vx.^2+md.initialization.vy.^2)/md.constants.yts;
ubx=md.initialization.vx/md.constants.yts;
uby=md.initialization.vy/md.constants.yts;

%compute basal drag (S.I.)
alpha2 = (N.^r).*(md.friction.coefficient.^2).*(ub.^(s-1));
b  =  alpha2.*ub;
bx = -alpha2.*ubx;
by = -alpha2.*uby;

%return magnitude of only one output is requested
if nargout==1
	bx = b;
end
