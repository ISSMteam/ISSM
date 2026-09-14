%Test Name: SquareSheetConstrainedFrictionEmulatorSSA2d
md=triangle(model(),'../Exp/Square.exp',200000.);
md=setmask(md,'','');
md=parameterize(md,'../Par/SquareSheetConstrained.par');
md=setflowequation(md,'SSA','all');

%use friction emulator
md.friction = frictionemulator();
md.friction.module_dir = '../Data/friction_emulator';
md.friction.pt_name = 'model_sqrt_weighted_standard.pt';
md.friction.py_name = 'friction_emulator.py';
md.friction.C = 20*ones(md.mesh.numberofvertices,1);

md.cluster=generic('name',oshostname(),'np',1);
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

