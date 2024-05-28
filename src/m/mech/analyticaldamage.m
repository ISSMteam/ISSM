function [damage,B,backstress]=analyticaldamage(md,varargin)
%ANALYTICALDAMAGE - compute damage for an ice shelf 
%
%	 This routine computes damage as a function of water/ice
%	 material properties, ice thickness, strain rate, and ice 
%	 rigidity.  The model must contain computed strain rates,
%	 either from observed or modeled ice velocities.
%
%   Available options:
%		- 'eq'			: analytical equation to use in the calculation.  Must be one of:
%								'Weertman1D' for a confined ice shelf free to flow in one direction
%								'Weertman2D' for an unconfined ice shelf free to spread in any direction
%								'Thomas' for a 2D ice shelf, taking into account full strain rate tensor (default)
%		- 'smoothing'	: the amount of smoothing to be applied to the strain rate data.
%								Type 'help averaging' for more information on its
%								usage. Defaults to 0.
%		- 'sigmab'		: a-priori backstress in opposition to the driving stress.
%								Defaults to 0 everywhere. 

%		- 'coordsys'	: coordinate system for calculating the strain rate
%							components. Must be one of: 
%				'longitudinal': x axis aligned along a flowline at every point (default)
%				'principal': x axis aligned along maximum principal strain rate
%					at every point
%				'xy': x and y axes same as in polar stereographic projection 
%
%   Return values:
%		'damage' which is truncated in the range [0,1-1e-9]
%
%	   'B' is the ice rigidity calculated assuming D=0 everywhere. 
%
%		'backstress' is the inferred backstress necessary to balance the
%		analytical solution (keeping damage within its appropriate limits, e.g. D
%		in [0,1]).
%
%   Usage:
%      [damage,B,backstress]=analyticaldamage(md,options)
%
%   Example:
%      [damage,B,backstress]=analyticaldamage(md,'eq','Weertman2D','smoothing',2,'sigmab',10e3,'coordsys','longitudinal');

% check inputs
if (nargin<1),
	help analyticaldamage
	error('bad usage');
end
if isempty(fieldnames(md.results)),
	error(['md.results.strainrate is not present.  Calculate using md=mechanicalproperties(md,vx,vy)']);
end
if dimension(md.mesh)~=2,
	error('only 2d model supported currently');
end
if any(md.flowequation.element_equation~=2),
	disp('Warning: the model has some non SSA elements. These will be treated like SSA elements');
end

% process options
options = pairoptions(varargin{:});
eq = getfieldvalue(options,'eq','Thomas');
smoothing = getfieldvalue(options,'smoothing',0);
sigmab = getfieldvalue(options,'sigmab',0);
coordsys = getfieldvalue(options,'coordsys','longitudinal');
if length(sigmab)==1,
	sigmab=sigmab*ones(md.mesh.numberofvertices,1);
end

[a,b,theta,ex]=thomasparams(md,'eq',eq,'smoothing',smoothing,'coordsys',coordsys);

% spreading stress
rhoi=md.materials.rho_ice;
rhow=md.materials.rho_water;
C=0.5*rhoi*md.constants.g*(1-rhoi/rhow);
T=C*md.geometry.thickness;

% rheology
B=md.materials.rheology_B;
n=averaging(md,md.materials.rheology_n,0);

D=1-(1+a+a.^2+b.^2).^((n-1)./(2*n))./abs(ex).^(1./n).*(T-sigmab)./B./(2+a)./sign(ex);

% D>1 where (2+a).*sign(ex)<0, compressive regions where high backstress needed
pos=find(D>1);
D(pos)=0;

backstress=zeros(md.mesh.numberofvertices,1);

% backstress to bring D down to one 
backstress(pos)=T(pos)-(1-D(pos)).*B(pos).*sign(ex(pos)).*(2+a(pos)).*abs(ex(pos)).^(1./n(pos))./(1+a(pos)+a(pos).^2).^((n(pos)-1)/2./n(pos));

pos=find(D<0);
mask=ismember(1:md.mesh.numberofvertices,pos);
D(pos)=0;

% backstress to bring negative damage to zero
backstress(pos)=T(pos)-(1-D(pos)).*B(pos).*sign(ex(pos)).*(2+a(pos)).*abs(ex(pos)).^(1./n(pos))./(1+a(pos)+a(pos).^2).^((n(pos)-1)/2./n(pos));

pos=find(backstress<0);
backstress(pos)=0;

% rigidity from Thomas relation for D=0 and backstress=0
B=sign(ex)./(2+a).*(1+a+a.^2).^((n-1)/2./n).*T./(abs(ex).^(1./n));
pos=find(B<0);
B(pos)=md.materials.rheology_B(pos);

damage=D;
