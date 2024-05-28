%Test Name: PigTherSteaSUPG
md=triangle(model(),'../Exp/Pig.exp',30000.);
md=setmask(md,'../Exp/PigShelves.exp','../Exp/PigIslands.exp');
md=parameterize(md,'../Par/Pig.par');
md=extrude(md,3,1.);
md=setflowequation(md,'HO','all');
md.thermal.stabilization=2;
md.cluster=generic('name',oshostname(),'np',3);
md.timestepping.time_step=0;
md.thermal.penalty_threshold=40;
md=solve(md,'Thermal');

%Fields and tolerances to track changes
field_names     ={'Temperature','BasalforcingsGroundediceMeltingRate'};
field_tolerances={1e-11,1e-11};
field_values={...
	(md.results.ThermalSolution.Temperature),...
	(md.results.ThermalSolution.BasalforcingsGroundediceMeltingRate),...
	};
