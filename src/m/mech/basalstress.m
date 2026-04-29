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
	%   See also: EFFECTIVEPRESSURE, PLOT_BASALDRAG

	if isprop(md.friction, 'coupling')
		coupling = md.friction.coupling;
	else
		warning(sprintf('md.friction.coupling is not found. Default coupling is set to 0.'));
		coupling = 0;
	end

	% calculate effective pressure using coupling definition in md.friction
	N = effectivepressure(md); % effective pressure (Pa)

	% compute sliding velocity
	ub  = sqrt(md.initialization.vx.^2+md.initialization.vy.^2)/md.constants.yts; % horizontal vel (m/s)
	ubx = md.initialization.vx/md.constants.yts; % vx (m/s)
	uby = md.initialization.vy/md.constants.yts; % vy (m/s)

	%compute basal drag (S.I.)
	switch(class(md.friction))
		case 'friction'
			%compute exponents
			s=averaging(md,1./md.friction.p,0);
			r=averaging(md,md.friction.q./md.friction.p,0);

			alpha2 = (N.^r).*(md.friction.coefficient.^2).*(ub.^(s-1));

		case 'frictionschoof'
			if any(N < 0)
				%NOTE: Negative values of effective pressure N return a complex number in alpha2. Treated here with minimum threshold.
				warning('Find effective pressure value N < 0. Enforcing minimum effective pressure of N_min = 0.1');
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
