md=model();

%Create mesh with roundmesh
md=roundmesh(md,750000,30000);

%Set mask
md=setmask(md,'','');

%Parameterize model
md=parameterize(md,'EISMINT.par');

%We extrude the model to have a 3d model
md=extrude(md,10,1.);

%Set ice flow approximation
md=setflowequation(md,'SIA','all');

%Create boundary conditions: zero velocity on the bed
pos=find(md.mesh.vertexonbase);
md.stressbalance.spcvx(pos)=0;
md.stressbalance.spcvy(pos)=0;
md.stressbalance.spcvz(pos)=0;

%Go Solve
md.cluster=generic('np',2);
md.verbose.convergence=1;
md=solve(md,'Stressbalance');
vel=DepthAverage(md,sqrt(md.results.StressbalanceSolution.Vx.^2+md.results.StressbalanceSolution.Vy.^2));

%Calculate analytical velocity
constant=0.3;
vx_obs=constant/2*md.mesh.x.*(md.geometry.thickness).^-1;
vy_obs=constant/2*md.mesh.y.*(md.geometry.thickness).^-1;
vel_obs=sqrt(vx_obs.^2+vy_obs.^2);
vel_obs=project2d(md,vel_obs,1);

plotmodel(md,...
	'data',vel    ,'view',2,'caxis',[0 200],'title','Modelled velocity',...
	'data',vel_obs,'view',2,'caxis',[0 200],'title','Analytical velocity',...
	'data',abs(vel-vel_obs)./(vel_obs+eps)*100,'caxis',[0 30],'view',2,'title','Relative misfit (%)');

subplot(2,2,4)
hold on;
plot(sqrt((md.mesh.x2d).^2+(md.mesh.y2d).^2),vel,'r.');
plot(sqrt((md.mesh.x2d).^2+(md.mesh.y2d).^2),vel_obs,'b.');
title('Analytical vs calculated velocity');
xlabel('distance to the center of the ice sheet (m)');
ylabel('velocity (m/yr)');
legend('calculated velocity','exact velocity');
axis([0 750000 0 200]);
hold off;
