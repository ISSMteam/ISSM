function [bx by b]=basalstress(md)
%BASALSTRESS - compute basal stress from basal drag and geometric information. 
%
%      Computes basal stress from geometric information and ice velocity in md.initialization.
%
%   Usage:
%      [bx by b]=basalstress(md);
%
%   See also: plot_basaldrag

%compute exponents
s=averaging(md,1./md.friction.p,0);
r=averaging(md,md.friction.q./md.friction.p,0);

%Compute effective pressure
switch(md.friction.coupling)
	case 0
		N = max(md.constants.g*(md.materials.rho_ice*md.geometry.thickness+md.materials.rho_water*md.geometry.base),0);
	case 3
		N = max(md.friction.effective_pressure, 0);
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
