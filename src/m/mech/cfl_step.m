function maxtime=cfl_step(md,vx,vy)
%CFL_STEP - return the maximum time step for the model in years
%
%   Dt < 0.5 / (u/Dx + v/Dy)
%
%   Usage:
%      maxtime=cfl_step(md,vx,vy);
%
%   Example:
%      dt=cfl_step(md,md.results.StressbalanceSolution.Vx,md.results.StressbalanceSolution.Vy)

%Check length of velocities 
if size(vx,1)~=md.mesh.numberofvertices & size(vy,1)~=md.mesh.numberofvertices,
	error('timesteps error message: size of velocity components must be the same as md.mesh.numberofvertices');
end

index=md.mesh.elements;
edgex=max(md.mesh.x(index),[],2)-min(md.mesh.x(index),[],2);
edgey=max(md.mesh.y(index),[],2)-min(md.mesh.y(index),[],2);
vx=max(abs(vx(index)),[],2);
vy=max(abs(vy(index)),[],2);

maxtime=1/2*min(1./(vx./edgex+vy./edgey));
