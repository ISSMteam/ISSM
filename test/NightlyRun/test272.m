%Test Name: SquareShelfCMDSSA2dDamage
md=triangle(model(),'../Exp/Square.exp',150000.);
md=setmask(md,'all','');
md.materials=matdamageice();
md=parameterize(md,'../Par/SquareShelf.par');
md.damage.isdamage=1;
md.damage.D=0.5*ones(md.mesh.numberofvertices,1);
md.damage.spcdamage=NaN*ones(md.mesh.numberofvertices,1);
md=setflowequation(md,'SSA','all');
md.verbose=verbose('control',true);

%control parameters
md.inversion.iscontrol=1;
md.inversion.control_parameters={'DamageDbar'};
md.inversion.min_parameters=zeros(md.mesh.numberofvertices,1);
md.inversion.max_parameters=(1-10^-13)*ones(md.mesh.numberofvertices,1);
md.inversion.nsteps=2;
md.inversion.cost_functions=101;
md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,1);
md.inversion.gradient_scaling=0.9*ones(md.inversion.nsteps,1);
md.inversion.maxiter_per_step=2.*ones(md.inversion.nsteps,1);
md.inversion.step_threshold=0.99*ones(md.inversion.nsteps,1);
md.inversion.vx_obs=md.initialization.vx; 
md.inversion.vy_obs=md.initialization.vy;

md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Stressbalance');

%Fields and tolerances to track changes
field_names     ={'Gradient','Misfits','DamageDbar','Pressure','Vel','Vx','Vy'};
field_tolerances={1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,1e-12};
field_values={...
   (md.results.StressbalanceSolution.Gradient1),...
   (md.results.StressbalanceSolution.J),...
   (md.results.StressbalanceSolution.DamageDbar),...
   (md.results.StressbalanceSolution.Pressure),...
   (md.results.StressbalanceSolution.Vel),...
   (md.results.StressbalanceSolution.Vx),...
   (md.results.StressbalanceSolution.Vy)
};
