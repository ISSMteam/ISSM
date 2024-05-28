%Test Name: 79NorthCMBalThic2dCG
md=triangle(model(),'../Exp/79North.exp',10000.);
md=setmask(md,'../Exp/79NorthShelf.exp','');
md=parameterize(md,'../Par/79North.par');
md=setflowequation(md,'SSA','all');

%control parameters
md.inversion.nsteps=2;
md.masstransport.stabilization=1;
md.inversion.iscontrol=1;
md.inversion.control_parameters={'BalancethicknessThickeningRate'};
md.inversion.thickness_obs=md.geometry.thickness;
md.inversion.min_parameters=-50.*ones(md.mesh.numberofvertices,1);
md.inversion.max_parameters=50.*ones(md.mesh.numberofvertices,1);
md.inversion.cost_functions=201;
md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,1);
md.inversion.gradient_scaling=10./md.constants.yts*ones(md.inversion.nsteps,1);
md.inversion.maxiter_per_step=4*ones(md.inversion.nsteps,1);
md.inversion.step_threshold=0.99*ones(md.inversion.nsteps,1);

md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Balancethickness');

%Fields and tolerances to track changes
field_names     ={'Gradient','Misfits','BalancethicknessThickeningRate','Thickness'};
field_tolerances={1e-12,1e-12,1e-12,1e-12,1e-12,1e-12};
field_values={...
	(md.results.BalancethicknessSolution.Gradient1),...
	md.results.BalancethicknessSolution.J,...
	(md.results.BalancethicknessSolution.BalancethicknessThickeningRate),...
	(md.results.BalancethicknessSolution.Thickness)
};
