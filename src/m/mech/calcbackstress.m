function backstress=calcbackstress(md,varargin)
%BACKSTRESSFROMINVERSION - compute ice shelf backstress  
%
%	 This routine computes backstress based on the analytical formalism of
%	 Thomas (1973) and Borstad et al. (2013) based on the ice rigidity,
%	 thickness, the densities of ice and seawater, and (optionally)
%	 damage. Strain rates must also be included, either from observed or
%	 modeled velocities.

%   Available options:
%		- 'smoothing'	: the amount of smoothing to be applied to the strain rate data.
%								Type 'help averaging' for more information on its
%								usage. Defaults to 0.
%		- 'coordsys'	: coordinate system for calculating the strain rate
%							components. Must be one of: 
%				'longitudinal': x axis aligned along a flowline at every point (default)
%				'principal': x axis aligned along maximum principal strain rate
%					at every point
%				'xy': x and y axes same as in polar stereographic projection 
%
%   Return values:
%		'backstress' is the inferred backstress based on the analytical
%		solution for ice shelf creep
%
%   Usage:
%      backstress=calcbackstress(md,options)
%
%   Example:
%      backstress=backstressfrominversion(md,'smoothing',2,'coordsys','longitudinal');

% check inputs
if (nargin<1),
	help backstressfrominversion
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
smoothing = getfieldvalue(options,'smoothing',0);
coordsys = getfieldvalue(options,'coordsys','longitudinal');
tempmask = getfieldvalue(options,'tempmask',false);

T=0.5*md.materials.rho_ice*md.constants.g*(1-md.materials.rho_ice/md.materials.rho_water)*md.geometry.thickness;
n=averaging(md,md.materials.rheology_n,0);
B=md.materials.rheology_B;
if md.damage.isdamage,
	D=md.damage.D
else
	D=0;
end

[a0,b0,theta0,ex0]=thomasparams(md,'eq','Thomas','smoothing',smoothing,'coordsys',coordsys);

% analytical backstress solution
%backstress=T-(1.-D).*B.*sign(ex0).*(2+a0).*abs(ex0).^(1./n)./((1+a0+a0.^2+b0.^2).^((n-1)/2./n));
backstress=1-((1.-D).*B.*sign(ex0).*(2+a0).*abs(ex0).^(1./n)./((1+a0+a0.^2+b0.^2).^((n-1)/2./n)))./T; %fix
backstress(find(backstress<0))=0;
