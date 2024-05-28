%Test Name: SquareSheetConstrainedStressSSA2d
md=triangle(model(),'../Exp/Square.exp',150000.);
md=setmask(md,'','');
md=parameterize(md,'../Par/SquareSheetConstrained.par');
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Stressbalance');

%Fields and tolerances to track changes
field_names     ={'Vx','Vy','Vel','Pressure'};
field_tolerances={1e-13,1e-13,1e-13,1e-13};
field_values={...
	(md.results.StressbalanceSolution.Vx),...
	(md.results.StressbalanceSolution.Vy),...
	(md.results.StressbalanceSolution.Vel),...
	(md.results.StressbalanceSolution.Pressure),...
	};
