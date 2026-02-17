function J=velocitymisfit(md)
%MISFIT - compute misfit
%
%   Usage:
%      J = velocitymisfit(md)
%
%   - modeled  velocities will be pulled from md.initialization.vx/vy
%   - observed velocities will be pulled from md.inversion.vx_obs/vy_obs


if dimension(md.mesh)==2,
	elements=md.mesh.elements;
	x=md.mesh.x;
	y=md.mesh.y;
	%vx=md.initialization.vx;
	%vy=md.initialization.vy;
	vx=md.results.StressbalanceSolution.Vx;
	vy=md.results.StressbalanceSolution.Vy;
	vx_obs=md.inversion.vx_obs;
	vy_obs=md.inversion.vy_obs;
else
	elements=md.mesh.elements2d;
	x=md.mesh.x2d;
	y=md.mesh.y2d;
	vx=project2d(md,md.initialization.vx,md.mesh.numberoflayers);
	vy=project2d(md,md.initialization.vy,md.mesh.numberoflayers);
	vx_obs=project2d(md,md.inversion.vx_obs,md.mesh.numberoflayers);
	vy_obs=project2d(md,md.inversion.vy_obs,md.mesh.numberoflayers);
end

%compute areas;
areas=GetAreas(elements,x,y);

%compute delta v on elements
deltav=1/2*(   (vx-vx_obs).^2+(vy-vy_obs).^2)/md.constants.yts^2;
deltav_elem=deltav(elements)*[1;1;1]/3;

%compute misfit
J=sum(deltav_elem.*areas);
