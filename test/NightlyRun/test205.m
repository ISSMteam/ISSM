%Test Name: SquareShelfStressMHOPenalties
md=triangle(model(),'../Exp/Square.exp',150000.);
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareShelf.par');
md=extrude(md,3,2.);
md=setflowequation(md,'HO','../Exp/SquareHalfRight.exp','fill','SSA','coupling','penalties');
md.settings.solver_residue_threshold = 1.e-4;
md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Stressbalance');

%Fields and tolerances to track changes
field_names     ={'Vx','Vy','Vz','Vel','Pressure'};
field_tolerances={2e-05,2e-05,1e-05,1e-05,1e-05};
field_values={...
	(md.results.StressbalanceSolution.Vx),...
	(md.results.StressbalanceSolution.Vy),...
	(md.results.StressbalanceSolution.Vz),...
	(md.results.StressbalanceSolution.Vel),...
	(md.results.StressbalanceSolution.Pressure),...
	};
