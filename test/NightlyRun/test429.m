%Test Name: SquareSheetShelfStressHONewton
md=triangle(model(),'../Exp/Square.exp',180000.);
md=setmask(md,'../Exp/SquareShelf.exp','');
md=parameterize(md,'../Par/SquareSheetShelf.par');
md=extrude(md,5,1.);
md=setflowequation(md,'HO','all');
md.stressbalance.isnewton=1;
md.stressbalance.restol=0.0001;
md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Stressbalance');

%Fields and tolerances to track changes
field_names     ={'Vx','Vy','Vz','Vel','Pressure'};
field_tolerances={2e-06,1e-06,1e-06,1e-06,1e-06};
field_values={...
	(md.results.StressbalanceSolution.Vx),...
	(md.results.StressbalanceSolution.Vy),...
	(md.results.StressbalanceSolution.Vz),...
	(md.results.StressbalanceSolution.Vel),...
	(md.results.StressbalanceSolution.Pressure),...
	};
