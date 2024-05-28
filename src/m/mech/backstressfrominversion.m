function backstress=backstressfrominversion(md,varargin)
%BACKSTRESSFROMINVERSION - compute ice shelf backstress from inversion results 
%
%	 This routine computes backstress based on the analytical formalism of
%	 Thomas (1973) and Borstad et al. (2013).  The model must contain inversion
%	 results for ice rigidity.  Strain rates must also be included, either from
%	 observed or modeled velocities.  Ice rigidity B is assumed to be
%	 parameterized by the ice temperature in md.materials.rheology_B.
%
%   Available options:
%		- 'tempmask'	: mask the inverted rigidity to be no more than
%							appropriate for the temperature of the ice?  
%							Boolean, defaults to false.
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
%      backstress=backstressfrominversion(md,options)
%
%   Example:
%      backstress=backstressfrominversion(md,'smoothing',2,'coordsys','longitudinal','tempmask',true);

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
Bi=md.results.StressbalanceSolution.MaterialsRheologyBbar;

[a0,b0,theta0,ex0]=thomasparams(md,'eq','Thomas','smoothing',smoothing,'coordsys',coordsys);

if tempmask
	Bi=md.results.StressbalanceSolution.MaterialsRheologyBbar;
   pos=find(Bi>md.materials.rheology_B);
   Bi(pos)=md.materials.rheology_B(pos);
end

% analytical backstress solution
backstress=T-Bi.*sign(ex0).*(2+a0).*abs(ex0).^(1./n)./((1+a0+a0.^2+b0.^2).^((n-1)/2./n));
backstress(find(backstress<0))=0;
