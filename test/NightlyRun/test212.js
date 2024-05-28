//Test Name: SquareShelfCMBSSA2d
var md = new model();
triangle(md,square[0],200000.);
setmask(md,'all','');
parameterize(md);
setflowequation(md,'SSA','all');

//control parameters
md.inversion.iscontrol=1;
md.inversion.control_parameters=['MaterialsRheologyBbar'];
md.inversion.min_parameters=Math.pow(10,6)*ones(md.mesh.numberofvertices,1);
md.inversion.max_parameters=2*Math.pow(10,9)*ones(md.mesh.numberofvertices,1);
md.inversion.nsteps=2;
md.inversion.cost_functions=101;
md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,1);
md.inversion.gradient_scaling=Math.pow(10,7)*ones(md.inversion.nsteps,1);
md.inversion.maxiter_per_step=2*ones(md.inversion.nsteps,1);
md.inversion.step_threshold=0.3*ones(md.inversion.nsteps,1);
md.inversion.vx_obs=md.initialization.vx; md.inversion.vy_obs=md.initialization.vy;
md.verbose.control=true;

//md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Stressbalance');

//Fields and tolerances to track changes
field_names     =['Gradient','Misfits','MaterialsRheologyBbar','Pressure','Vel','Vx','Vy'];
field_tolerances=[1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13];
field_values=[
	(md.results.StressbalanceSolution[0].Gradient1),
	(md.results.StressbalanceSolution[0].J),
	(md.results.StressbalanceSolution[0].MaterialsRheologyBbar),
	(md.results.StressbalanceSolution[0].Pressure),
	(md.results.StressbalanceSolution[0].Vel),
	(md.results.StressbalanceSolution[0].Vx),
	(md.results.StressbalanceSolution[0].Vy)
];
