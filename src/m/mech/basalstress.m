function [bx by b]=basalstress(md)
%BASALSTRESS - compute basal stress from basal drag and geometric information. 
%
%      Computes basal stress from geometric information and ice velocity in md.initialization. Follows the basal stress definition in "src/c/classes/Loads/Friction.cpp", lines 1102-1136.
%
%   Usage:
%      [bx by b]=basalstress(md);
%
%   See also: plot_basaldrag

%Compute effective pressure
g         = md.constants.g;
rho_ice   = md.materials.rho_ice;
rho_water = md.materials.rho_water;

sealevel = 0;
p_ice = g*rho_ice*md.geometry.thickness;

if isprop(md.friction, 'coupling')
	coupling = md.friction.coupling;
else
	warning(sprintf('md.friction.coupling is not found. Default coupling is set to 0.'));
	coupling = 0;
end

switch(coupling)
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
switch(class(md.friction))
	case 'friction'
		%compute exponents
		s=averaging(md,1./md.friction.p,0);
		r=averaging(md,md.friction.q./md.friction.p,0);

		alpha2 = (N.^r).*(md.friction.coefficient.^2).*(ub.^(s-1));

	case 'frictionschoof'
		if any(N < 0)
			%NOTE: Sometimes, negative value of effective pressure N gives image number in alpha2. To prevent the image value in alpha2, we use small values.
			warning('Find effective pressure value < 0. Treating minimum effective value with 0.1');
			N = max(N, 0.1);
		end
		m=averaging(md,md.friction.m,0);
		C=averaging(md,md.friction.C,0);
		Cmax=averaging(md,md.friction.Cmax,0);
		
		alpha2 = (C.^2 .* ub.^(m-1))./(1 + (C.^2./(Cmax.*N)).^(1./m).*ub).^(m);

	case 'frictionweertman'
		m = averaging(md,md.friction.m,0);
		C = md.friction.C;
		alpha2 = C.^2 .* ub.^(1./m-1);

	otherwise
		error('friction class not supported yet');
end
b  =  alpha2.*ub;
bx = -alpha2.*ubx;
by = -alpha2.*uby;

%return magnitude of only one output is requested
if nargout==1
	bx = b;
end
