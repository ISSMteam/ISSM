%Test Name: PigCMDragHO
md=triangle(model(),'../Exp/Pig.exp',20000.);
md=setmask(md,'../Exp/PigShelves.exp','../Exp/PigIslands.exp');
md=parameterize(md,'../Par/Pig.par');
md=extrude(md,4,1.);
md=setflowequation(md,'HO','all');

%control parameters
md.inversion.iscontrol=1;
md.inversion.control_parameters={'FrictionCoefficient'};
md.inversion.min_parameters=1.*ones(md.mesh.numberofvertices,1);
md.inversion.max_parameters=200.*ones(md.mesh.numberofvertices,1);
md.inversion.nsteps=2;
md.inversion.cost_functions=[103  501];
md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,2); md.inversion.cost_functions_coefficients(:,2)=2.*10^-7;
md.inversion.gradient_scaling=3.*ones(md.inversion.nsteps,1);
md.inversion.maxiter_per_step=2*ones(md.inversion.nsteps,1);
md.inversion.step_threshold=0.99*ones(md.inversion.nsteps,1);
md.inversion.vx_obs=md.initialization.vx; md.inversion.vy_obs=md.initialization.vy;

md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Stressbalance');

%Fields and tolerances to track changes
field_names     ={'Gradient','Misfits','FrictionCoefficient','Pressure','Vel','Vx','Vy'};
field_tolerances={1e-11,1e-11,1e-11,1e-11,1e-11,1e-11,1e-11,1e-11,1e-11,1e-11};
field_values={...
	(md.results.StressbalanceSolution.Gradient1),...
	md.results.StressbalanceSolution.J,...
	(md.results.StressbalanceSolution.FrictionCoefficient),...
	(md.results.StressbalanceSolution.Pressure),...
	(md.results.StressbalanceSolution.Vel),...
	(md.results.StressbalanceSolution.Vx),...
	(md.results.StressbalanceSolution.Vy)
};
