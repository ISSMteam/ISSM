%Test Name: SquareSheetShelfStressSSAFS3dTiling
md=triangle(model(),'../Exp/Square.exp',180000.);
md=setmask(md,'../Exp/SquareShelf.exp','');
md=parameterize(md,'../Par/SquareSheetShelf.par');
md=extrude(md,5,1.);
md=setflowequation(md,'FS','../Exp/SquareHalfRight.exp','fill','SSA');
md.cluster=generic('name',oshostname(),'np',3);
md.stressbalance.reltol=0.4;
md=solve(md,'Stressbalance');

%Fields and tolerances to track changes
field_names     ={'Vx','Vy','Vz','Vel','Pressure'};
field_tolerances={4e-07,4e-07,2e-06,4e-07,5e-07};
field_values={...
	(md.results.StressbalanceSolution.Vx),...
	(md.results.StressbalanceSolution.Vy),...
	(md.results.StressbalanceSolution.Vz),...
	(md.results.StressbalanceSolution.Vel),...
	(md.results.StressbalanceSolution.Pressure),...
	};
