%Test Name: PigCMBFSm1qn3
md=triangle(model(),'../Exp/Pig.exp',11000.);
md=setmask(md,'../Exp/PigShelves.exp','../Exp/PigIslands.exp');
md=parameterize(md,'../Par/Pig.par');

%impose hydrostatic equilibrium (required by Stokes)
md.geometry.base=-md.materials.rho_ice/md.materials.rho_water*md.geometry.thickness;
md.geometry.surface=md.geometry.base+md.geometry.thickness;
md=extrude(md,3,1.);
md=setflowequation(md,'FS','all');
md=extract(md,md.mask.ocean_levelset<0.);

%control parameters
md.inversion.iscontrol=1;
md.inversion.control_parameters={'MaterialsRheologyBbar'};
md.inversion.min_parameters=10.^6*ones(md.mesh.numberofvertices,1);
md.inversion.max_parameters=2.*10^9*ones(md.mesh.numberofvertices,1);
md.inversion.nsteps=2;
md.inversion.cost_functions=[101 502];
md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,2);
md.inversion.gradient_scaling=10.^8*ones(md.inversion.nsteps,1);
md.inversion.maxiter_per_step=2.*ones(md.inversion.nsteps,1);
md.inversion.step_threshold=0.99*ones(md.inversion.nsteps,1);
md.inversion.vx_obs=md.initialization.vx; md.inversion.vy_obs=md.initialization.vy;

md.inversion=m1qn3inversion(md.inversion);

md.cluster=generic('name',oshostname(),'np',1);
md=solve(md,'Stressbalance');

%Fields and tolerances to track changes
field_names     ={'Gradient','Misfits','MaterialsRheologyB','Pressure','Vel','Vx','Vy'};
field_tolerances={6e-11,5e-11,5e-10,1e-09,2e-11,5e-11,2e-11};
field_values={...
	(md.results.StressbalanceSolution.Gradient1),...
	(md.results.StressbalanceSolution.J),...
	(md.results.StressbalanceSolution.MaterialsRheologyBbar),...
	(md.results.StressbalanceSolution.Pressure),...
	(md.results.StressbalanceSolution.Vel),...
	(md.results.StressbalanceSolution.Vx),...
	(md.results.StressbalanceSolution.Vy)
};
