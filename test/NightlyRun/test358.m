%Test Name: SquareSheetConstrainedCMDragRegCoulombSSA2d
md=triangle(model(),'../Exp/Square.exp',200000.);
md=setmask(md,'','');
md=parameterize(md,'../Par/SquareSheetConstrained.par');
md=setflowequation(md,'SSA','all');

%use Regularized Coulomb's law
md.friction = frictionregcoulomb();
md.friction.m = 3.0*ones(md.mesh.numberofelements,1);
md.friction.u0 = 2000; %m/yr
md.friction.C = 200*ones(md.mesh.numberofvertices,1);
	
%control parameters
md.inversion.iscontrol=1;
md.inversion.control_parameters={'FrictionC'};
md.inversion.min_parameters=1.*ones(md.mesh.numberofvertices,1);
md.inversion.max_parameters=10000.*ones(md.mesh.numberofvertices,1);
md.inversion.nsteps=2;
md.inversion.cost_functions=[102  501];
md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,2); md.inversion.cost_functions_coefficients(:,2)=2.*10^-7;
md.inversion.gradient_scaling=3.*ones(md.inversion.nsteps,1);
md.inversion.maxiter_per_step=2*ones(md.inversion.nsteps,1);
md.inversion.step_threshold=0.3*ones(md.inversion.nsteps,1);
md.inversion.vx_obs=md.initialization.vx; md.inversion.vy_obs=md.initialization.vy;

md.verbose = verbose('all');
md.debug.valgrind = 0;

md.cluster=generic('name',oshostname(),'np',1);
md=solve(md,'Stressbalance');

%Fields and tolerances to track changes
field_names     ={'Gradient','Misfits','FrictionC','Pressure','Vel','Vx','Vy'};
field_tolerances={1e-12,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13};
field_values={...
	(md.results.StressbalanceSolution.Gradient1),...
	(md.results.StressbalanceSolution.J),...
	(md.results.StressbalanceSolution.FrictionC),...
	(md.results.StressbalanceSolution.Pressure),...
	(md.results.StressbalanceSolution.Vel),...
	(md.results.StressbalanceSolution.Vx),...
	(md.results.StressbalanceSolution.Vy)
};
