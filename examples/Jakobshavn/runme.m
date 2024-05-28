steps=[1];

if any(steps==1) 
	disp('   Step 1: Mesh creation');
	md=triangle(model,'Domain.exp',2000);

	%Get observed velocity field on mesh nodes
	ncdata='../Data/Greenland_5km_dev1.2.nc';
	if ~exist(ncdata,'file'),
		error('File Greenland_5km_dev1.2.nc not downloaded in Data Directory');
	end
	x1		= ncread(ncdata,'x1');
	y1		= ncread(ncdata,'y1');
	velx	= ncread(ncdata,'surfvelx');
	vely	= ncread(ncdata,'surfvely');
	vx		= InterpFromGridToMesh(x1,y1,velx',md.mesh.x,md.mesh.y,0);
	vy		= InterpFromGridToMesh(x1,y1,vely',md.mesh.x,md.mesh.y,0);
	vel		= sqrt(vx.^2+vy.^2);

	%refine mesh using surface velocities as metric
	md=bamg(md,'hmin',1200,'hmax',15000,'field',vel,'err',5);
	[md.mesh.lat,md.mesh.long]  = xy2ll(md.mesh.x,md.mesh.y,+1,39,71);

	save JksMesh md
end 

if any(steps==2) 
	disp('   Step 2: Parameterization');
	md=loadmodel('JksMesh');

	md=setmask(md,'','');
	md=parameterize(md,'Jks.par');

	save JksPar md
end 

if any(steps==3) 
	disp('   Step 3: Control method friction');
	md=loadmodel('JksPar');

	md=setflowequation(md,'SSA','all');

	%Control general
	md.inversion.iscontrol=1;
	md.inversion.nsteps=20;
	md.inversion.step_threshold=0.99*ones(md.inversion.nsteps,1);
	md.inversion.maxiter_per_step=5*ones(md.inversion.nsteps,1);
	md.verbose=verbose('solution',true,'control',true);

	%Cost functions
	md.inversion.cost_functions=[101 103];
	md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,2);
	md.inversion.cost_functions_coefficients(:,1)=40;
	md.inversion.cost_functions_coefficients(:,2)=1;

	%Controls
	md.inversion.control_parameters={'FrictionCoefficient'};
	md.inversion.gradient_scaling(1:md.inversion.nsteps)=30;
	md.inversion.min_parameters=1*ones(md.mesh.numberofvertices,1);
	md.inversion.max_parameters=200*ones(md.mesh.numberofvertices,1);

	%Additional parameters
	md.stressbalance.restol=0.01;
	md.stressbalance.reltol=0.1;
	md.stressbalance.abstol=NaN;

	%Go solve
	md.cluster=generic('name',oshostname,'np',4);
	md=solve(md,'Stressbalance');

	save JksControl md
end 

if any(steps==4) 
	disp('   Plotting')
	md=loadmodel('JksControl');

	plotmodel(md,'unit#all','km','axis#all','equal',...
		'data',md.inversion.vel_obs,'title','Observed velocity',...
		'data',md.results.StressbalanceSolution.Vel,'title','Modeled Velocity',...
		'colorbar#1','off','colorbartitle#2','(m/yr)',...
		'caxis#1-2',[0,7000],...
		'data',md.geometry.base,'title','Base elevation',...
		'data',md.results.StressbalanceSolution.FrictionCoefficient,...
		'title','Friction Coefficient',...
		'colorbartitle#3','(m)');
end 
