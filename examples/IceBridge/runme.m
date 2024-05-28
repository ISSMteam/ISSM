clear all;
steps=[1];

%Location of SeaRISE dataset
ncdata='../Data/Greenland_5km_dev1.2.nc';

if any(steps==1) 
	disp('   Step 1: Mesh creation');

	%Generate initial uniform mesh (resolution = 20000 m)
	md=triangle(model,'./DomainOutline.exp',20000);

	% Get velocities (Note: You can use ncdisp('file') to see an ncdump)
	x1   = ncread(ncdata,'x1');
	y1   = ncread(ncdata,'y1');
	velx = ncread(ncdata,'surfvelx');
	vely = ncread(ncdata,'surfvely');
	vx   = InterpFromGridToMesh(x1,y1,velx',md.mesh.x,md.mesh.y,0);
	vy   = InterpFromGridToMesh(x1,y1,vely',md.mesh.x,md.mesh.y,0);
	vel  = sqrt(vx.^2+vy.^2);

	%Mesh greenland without refinement in Jak basin
	md=bamg(md,'hmax',400000,'hmin',5000,'gradation',1.7,'field',vel,'err',8);
	return;

	%Refine mesh in the region of Jakobshavn (resolution = 3000 m)
	hmaxVertices=NaN*ones(md.mesh.numberofvertices,1);
	in=ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,...
		'Jak_outline.exp','node',1);
	hmaxVertices(find(in))=3000;
	md=bamg(md,'hmax',400000,'hmin',5000,'gradation',1.7,'field',vel,...
		'err',8,'hmaxVertices',hmaxVertices);

	%convert x,y coordinates (Polar stereo) to lat/lon
	[md.mesh.lat,md.mesh.long]=xy2ll(md.mesh.x,md.mesh.y,+1,39,71);

	save ./Models/Greenland.Mesh_generation md;
end 

if any(steps==2) 
	disp('   Step 2: Parameterization');
	md = loadmodel('./Models/Greenland.Mesh_generation');

	md = setmask(md,'','');
	md = parameterize(md,'./Greenland.par');
	md = setflowequation(md,'SSA','all');

	save ./Models/Greenland.Parameterization2 md; 
end 

if any(steps==3) 
	disp('   Step 3: Control method friction');
	md = loadmodel('./Models/Greenland.Parameterization2');

	%Control general
	md.inversion.iscontrol=1;
	md.inversion.nsteps=30;
	md.inversion.step_threshold=0.99*ones(md.inversion.nsteps,1);
	md.inversion.maxiter_per_step=5*ones(md.inversion.nsteps,1);
	md.verbose=verbose('solution',true,'control',true);

	%Cost functions
	md.inversion.cost_functions=[101 103 501];
	md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,3);
	md.inversion.cost_functions_coefficients(:,1)=350;
	md.inversion.cost_functions_coefficients(:,2)=0.6;
	md.inversion.cost_functions_coefficients(:,3)=2e-6;

	%Controls
	md.inversion.control_parameters={'FrictionCoefficient'};
	md.inversion.gradient_scaling(1:md.inversion.nsteps)=50;
	md.inversion.min_parameters=1*ones(md.mesh.numberofvertices,1);
	md.inversion.max_parameters=200*ones(md.mesh.numberofvertices,1);
	in=ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,...
		'data_gaps.exp','node',1);
	md.inversion.cost_functions_coefficients(find(in),1)=0.0;
	md.inversion.cost_functions_coefficients(find(in),2)=0.0;

	%Additional parameters
	md.stressbalance.restol=0.01;
	md.stressbalance.reltol=0.1;
	md.stressbalance.abstol=NaN;

	%Go solve
	md.cluster=generic('name',oshostname,'np',2);

	md.verbose=verbose('solution',true,'control',true);
	md=solve(md,'Stressbalance');

	%Update model friction fields accordingly
	md.friction.coefficient=md.results.StressbalanceSolution.FrictionCoefficient;

	save ./Models/Greenland.Control_drag md; 
end 

if any(steps==4) 
	disp('   Step 4: Transient run');
	md = loadmodel('./Models/Greenland.Control_drag');

	%Set surface mass balance
	x1  = ncread(ncdata,'x1');
	y1  = ncread(ncdata,'y1');
	smb = ncread(ncdata,'smb');
	smb = InterpFromGridToMesh(x1,y1,smb',md.mesh.x,md.mesh.y,0)*1000/md.materials.rho_ice;
	smb = [smb smb smb-1.0];
	md.smb.mass_balance = [smb;1 10 20];

	%Set transient options, run for 20 years, saving every year
	md.timestepping.time_step=0.2; %This must be reduced for finer resolutions
	md.timestepping.final_time=20;
	md.settings.output_frequency=1;

	%Additional options
	md.inversion.iscontrol=0;
	md.transient.requested_outputs={'IceVolume','TotalSmb','SmbMassBalance'};
	md.verbose=verbose('solution',true,'module',true,'convergence',true);

	%Go solve
	md.cluster=generic('name',oshostname,'np',2);
	md=solve(md,'Transient');

	save ./Models/Greenland.Transient md; 
end 

if any(steps==5) 
	disp('   Step 5: Plotting'); 
	md = loadmodel('./Models/Greenland.Transient');

	%Planview plots
	plotmodel(md,'data',md.results.TransientSolution(end).Vel,'caxis',[1e-1 6000],...
		'log', 10, 'title', 'Velocity (m/y)','gridded',1, ...
		'data', md.results.TransientSolution(end).SmbMassBalance, ...
		'title', 'Surface mass balance (m/y)','gridded',1, ...
		'data',md.results.TransientSolution(end).Thickness,...
		'title','Thickness (m)','gridded',1, ...
		'data',md.results.TransientSolution(end).Surface, ...
		'title', 'Surface (m)','gridded',1);

	%Line Plots
	figure
	time_plot=md.results.TransientSolution(1).time:md.timestepping.time_step:md.timestepping.final_time;

	%Plot surface mass balance
	surfmb=[]; for i=1:100; surfmb=[surfmb ...
		md.results.TransientSolution(i).SmbMassBalance]; end
	subplot(3,1,1); plot(time_plot,mean(surfmb)); title('Mean Surface mass balance');

	%Plot velocity
	vel=[]; for i=1:100; vel=[vel md.results.TransientSolution(i).Vel]; end
	subplot(3,1,2); plot(time_plot,mean(vel)); title('Mean Velocity');

	%Plot Volume
	volume=[]; for i=1:100; volume=[volume md.results.TransientSolution(i).IceVolume]; end
	subplot(3,1,3); plot(time_plot,volume); title('Ice Volume');
	xlabel('years')
end 
