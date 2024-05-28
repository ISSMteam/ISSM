%Test Name: SquareShelfConstrainedStressFSAdolc
md=triangle(model(),'../Exp/Square.exp',180000.);
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareShelfConstrained.par');
md=extrude(md,3,2.);
md=setflowequation(md,'FS','all');
md.cluster=generic('name',oshostname(),'np',1);
md.autodiff.isautodiff=true;
md.toolkits.DefaultAnalysis=issmgslsolver();
md.verbose=verbose('autodiff',true);
md=solve(md,'Stressbalance');

%Fields and tolerances to track changes
field_names     ={'Vx','Vy','Vz','Vel','Pressure'};
field_tolerances={1e-08,1e-08,3e-07,1e-08,1e-08};
field_values={...
	(md.results.StressbalanceSolution.Vx),...
	(md.results.StressbalanceSolution.Vy),...
	(md.results.StressbalanceSolution.Vz),...
	(md.results.StressbalanceSolution.Vel),...
	(md.results.StressbalanceSolution.Pressure),...
	};
