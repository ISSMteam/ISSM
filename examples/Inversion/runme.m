steps=[1];

if any(steps==1) 
	%Generate observations
	md = model;
	md = triangle(md,'DomainOutline.exp',100000);
	md = setmask(md,'all','');
	md = parameterize(md,'Square.par');
	md = setflowequation(md,'SSA','all');
	md.cluster = generic('np',2);
	md = solve(md,'Stressbalance');
	plotmodel(md,'axis#all','tight','data',md.materials.rheology_B,'caxis',[1.3 1.9]*10^8,'title','"True" B',...
		'data',md.results.StressbalanceSolution.Vel,'title','"observed velocities"')
	save model1 md
end 

if any(steps==2) 
	%Modify rheology, now constant
	loadmodel('model1.mat');
	md.materials.rheology_B(:) = 1.8*10^8;

	%results of previous run are taken as observations
	md.inversion=m1qn3inversion();
	md.inversion.vx_obs		= md.results.StressbalanceSolution.Vx;
	md.inversion.vy_obs		= md.results.StressbalanceSolution.Vy;
	md.inversion.vel_obs	= md.results.StressbalanceSolution.Vel;

	md = solve(md,'Stressbalance');
	plotmodel(md,'axis#all','tight','data',md.materials.rheology_B,'caxis',[1.3 1.9]*10^8,'title','B first guess',...
		'data',md.results.StressbalanceSolution.Vel,'title','modeled velocities')
	save model2 md
end 

if any(steps==3) 
	%invert for ice rigidity
	loadmodel('model2.mat');

	%Set up inversion parameters
	maxsteps = 20;
	md.inversion.iscontrol = 1;
	md.inversion.control_parameters = {'MaterialsRheologyBbar'};
	md.inversion.maxsteps = maxsteps;
	md.inversion.cost_functions = 101;
	md.inversion.cost_functions_coefficients = ones(md.mesh.numberofvertices,1);
	md.inversion.min_parameters = cuffey(273)*ones(md.mesh.numberofvertices,1);
	md.inversion.max_parameters = cuffey(200)*ones(md.mesh.numberofvertices,1);

	%Go solve!
	md.verbose=verbose(0);
	md=solve(md,'Stressbalance');
	plotmodel(md,'axis#all','tight','data',md.results.StressbalanceSolution.MaterialsRheologyBbar,'caxis',[1.3 1.9]*10^8,'title','inferred B',...
		'data',md.results.StressbalanceSolution.Vel,'title','modeled velocities')
end 

if any(steps==4) 
	%invert for ice rigidity
	loadmodel('model2.mat');

	%Set up inversion parameters
	maxsteps = 20;
	md.inversion.iscontrol = 1;
	md.inversion.control_parameters = {'MaterialsRheologyBbar'};
	md.inversion.maxsteps = maxsteps;
	md.inversion.cost_functions = [101 502];
	md.inversion.cost_functions_coefficients		= ones(md.mesh.numberofvertices,1);
	md.inversion.cost_functions_coefficients(:,2)	= 10^-16*ones(md.mesh.numberofvertices,1);
	md.inversion.min_parameters = cuffey(273)*ones(md.mesh.numberofvertices,1);
	md.inversion.max_parameters = cuffey(200)*ones(md.mesh.numberofvertices,1);

	%Go solve!
	md.verbose=verbose(0);
	md=solve(md,'Stressbalance');
	plotmodel(md,'axis#all','tight','data',md.results.StressbalanceSolution.MaterialsRheologyBbar,'caxis',[ 1.3 1.9]*10^8,'title','inferred B',...
		'data',md.results.StressbalanceSolution.Vel,'title','modeled velocities')
end 
