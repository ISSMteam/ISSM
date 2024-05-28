function damage=damagefrominversion(md)
%DAMAGEFROMINVERSION - compute ice shelf damage from inversion results
%
%	This routine computes damage based on the analytical formalism of Borstad et
%	al. (2013).  The model must contain inversion results for ice rigidity.  Ice
%	rigidity B is assumed to be parameterized by the ice temperature in
%	md.materials.rheology_B. 
%
%	Usage:
%		damage=damagefrominversion(md)
%
%	Example:
%		damage=damagefrominversion(md)

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

damage=zeros(md.mesh.numberofvertices,1);
% Damage where Bi softer than B(T)
pos=find(md.results.StressbalanceSolution.MaterialsRheologyBbar<md.materials.rheology_B);
damage(pos)=1-md.results.StressbalanceSolution.MaterialsRheologyBbar(pos)./md.materials.rheology_B(pos);
