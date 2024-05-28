%Test Name: SquareSheetShelfCMDragSteaSSA3d
md=triangle(model(),'../Exp/Square.exp',170000.);
md=setmask(md,'../Exp/SquareShelf.exp','');
md=parameterize(md,'../Par/SquareSheetShelf.par');
md=extrude(md,3,1.);
md=setflowequation(md,'SSA','all');

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
md.inversion.step_threshold=0.3*ones(md.inversion.nsteps,1);
md.timestepping.time_step=0.;
md.inversion.vx_obs=md.initialization.vx; md.inversion.vy_obs=md.initialization.vy;

md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Steadystate');

%Fields and tolerances to track changes
field_names     ={'Gradient','Misfits','FrictionCoefficient','Pressure','Vel','Vx','Vy','Vz','Temperature','BasalforcingsGroundediceMeltingRate'};
field_tolerances={1e-09,1e-9,3e-10,1e-13,1e-09,1e-09,1e-09,1e-8,1e-09,1e-6};
field_values={...
	(md.results.SteadystateSolution.Gradient1),...
	md.results.SteadystateSolution.J,...
	(md.results.SteadystateSolution.FrictionCoefficient),...
	(md.results.SteadystateSolution.Pressure),...
	(md.results.SteadystateSolution.Vel),...
	(md.results.SteadystateSolution.Vx),...
	(md.results.SteadystateSolution.Vy),...
	(md.results.SteadystateSolution.Vz),...
	(md.results.SteadystateSolution.Temperature),...
	(md.results.SteadystateSolution.BasalforcingsGroundediceMeltingRate)
};
