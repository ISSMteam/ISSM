%Test Name: PigBalVel2
md=triangle(model(),'../Exp/Pig.exp',20000.);
md=setmask(md,'../Exp/PigShelves.exp','../Exp/PigIslands.exp');
md=parameterize(md,'../Par/Pig.par');
md.initialization.vx(:)=0;
md.initialization.vy(:)=0;
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Balancevelocity');

%Fields and tolerances to track changes
field_names     ={'DrivingStressX','DrivingStressX','Vel'};
field_tolerances={1e-13,1e-13,1e-13};
field_values={...
	(md.results.BalancevelocitySolution.DrivingStressX),...
	(md.results.BalancevelocitySolution.DrivingStressY),...
	(md.results.BalancevelocitySolution.Vel),...
	};
