%Test Name: SquareShelfStressSSA2dRift
md=triangle(model(),'../Exp/SquareHole.exp','../Exp/Rifts.exp',50000.);
md=meshprocessrifts(md,'../Exp/Square.exp');
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareShelf2.par');
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',3);

%rift settings
md.rifts.riftstruct.fill='Melange';
md.rifts.riftstruct.fraction=0;
md.stressbalance.rift_penalty_lock=2;
md.stressbalance.rift_penalty_threshold=0;
md.rifts.riftstruct.fractionincrement=.1;
md=solve(md,'Stressbalance');

%Fields and tolerances to track changes
field_names     ={'Vx','Vy','Vel','Pressure'};
field_tolerances={9e-7,7e-8,9e-8,2e-11};
field_values={...
	(md.results.StressbalanceSolution.Vx),...
	(md.results.StressbalanceSolution.Vy),...
	(md.results.StressbalanceSolution.Vel),...
	(md.results.StressbalanceSolution.Pressure),...
	};
