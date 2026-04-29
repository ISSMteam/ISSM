function [bx by b]=basalstress(md)
% BASALSTRESS - compute basal stress from friction law and geometric information. 
% Computes basal stress from basal sliding parametrization in md.friction and geometry and ice velocity in md.initialization. Follows the basal stress definition in "src/c/classes/Loads/Friction.cpp", lines 1102-1136. 
%
% USEAGE:
%      b         = basalstress(md); % one argout returns the scalar magnitude				
%      [bx by b] = basalstress(md); % multiple argout returns the horizontal vector components
% INPUT:
%   md		ISSM model from which to take md.friction and md.initialization.
% OUTPUT:
%   bx		x component of basal stress
%   by		y component of basal stress
%   b			scalar magnitude of basal stress
%
%   See also: plot_basaldrag

%Check md.friction class
if ~(isa(md.friction,'friction') | isa(md.friction,'frictionschoof'))
	error('Error: md.friction only supports "friction.m", "frictionschoof.m" classes.');
end

% calculate effective pressure using coupling definition in md.friction (Pa)
sealevel= 0; % reference sea level (m)
p_ice   = md.constants.g * md.materials.rho_ice * md.geometry.thickness; % ice pressure (Pa)
p_water = md.constants.g * md.materials.rho_water * (sealevel-md.geometry.base); % water pressure (Pa)
switch(md.friction.coupling)
   case 0 % uniform sheet (negative water pressure ok, default)
      N = p_ice-p_water;
   case 1 % effective pressure is equal to the overburden pressure
      N = p_ice;
   case 2 % uniform sheet (water pressure >= 0)
      p_water=max(p_water,0.0);
      N = p_ice-p_water;
   case 3 % Use effective pressure prescrived in md.friction.effective_pressure
      N = max(md.friction.effective_pressure, 0);
   case 4 % Use effective pressure dynamically calculated by the hydrology model (i.e., fully coupled)
      error('md.friction.coupling=4 is not supported yet.');
   otherwise
      error('not supported yet');
end

% compute sliding velocity
ub  = sqrt(md.initialization.vx.^2+md.initialization.vy.^2)/md.constants.yts; % horizontal vel (m/s)
ubx = md.initialization.vx/md.constants.yts; % vx (m/s)
uby = md.initialization.vy/md.constants.yts; % vy (m/s)

%compute basal drag (S.I.)
switch(class(md.friction))
	case 'friction'
		% compute on the vertices if on the elements
		s=averaging(md,1./md.friction.p,0);
		r=averaging(md,md.friction.q./md.friction.p,0);

		alpha2 = (N.^r).*(md.friction.coefficient.^2).*(ub.^(s-1));
	case 'frictionschoof'
		if any(N < 0)
			%NOTE: Negative values of effective pressure N return a complex number in alpha2. Treated here with minimum threshold.
			warning('Find effective pressure value N < 0. Enforcing minimum effective pressure of N_min = 0.1');
			N = max(N, 0.1);
		end
		% compute on the vertices if on the elements
		m=averaging(md,md.friction.m,0);
		C=averaging(md,md.friction.C,0);
		Cmax=averaging(md,md.friction.Cmax,0);
		
		alpha2 = (C.^2 .* ub.^(m-1))./(1 + (C.^2./(Cmax.*N)).^(1./m).*ub).^(m);
	otherwise
		error('not supported yet');
end
b  =  alpha2.*ub;
bx = -alpha2.*ubx;
by = -alpha2.*uby;

%return magnitude of only one output is requested
if nargout==1
	bx = b;
end
