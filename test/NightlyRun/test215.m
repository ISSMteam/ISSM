%Test Name: SquareShelfCMBFS
md=triangle(model(),'../Exp/Square.exp',200000.);
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareShelf.par');
md=extrude(md,3,1.);
md=setflowequation(md,'FS','all');
md.settings.solver_residue_threshold = 1.e-4;

%control parameters
md.inversion.iscontrol=1;
md.inversion.control_parameters={'MaterialsRheologyBbar'};
md.inversion.min_parameters=10^6*ones(md.mesh.numberofvertices,1);
md.inversion.max_parameters=2*10^9*ones(md.mesh.numberofvertices,1);
md.inversion.nsteps=2;
md.inversion.cost_functions=101;
md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,1);
md.inversion.gradient_scaling=10^7*ones(md.inversion.nsteps,1);
md.inversion.maxiter_per_step=2*ones(md.inversion.nsteps,1);
md.inversion.step_threshold=0.3*ones(md.inversion.nsteps,1);
md.inversion.vx_obs=md.initialization.vx; md.inversion.vy_obs=md.initialization.vy;

md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Stressbalance');

%Fields and tolerances to track changes
field_names     ={'Gradient','Misfits','MaterialsRheologyBbar','Pressure','Vel','Vx','Vy'};
field_tolerances={4.6e-08,1e-08,2e-08,2e-09,3e-09,2e-08,2e-08};
field_values={...
	(md.results.StressbalanceSolution.Gradient1),...
	(md.results.StressbalanceSolution.J),...
	(md.results.StressbalanceSolution.MaterialsRheologyBbar),...
	(md.results.StressbalanceSolution.Pressure),...
	(md.results.StressbalanceSolution.Vel),...
	(md.results.StressbalanceSolution.Vx),...
	(md.results.StressbalanceSolution.Vy)
};
