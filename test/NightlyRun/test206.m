%Test Name: SquareShelfTherStea
md=triangle(model(),'../Exp/Square.exp',180000.);
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareShelf.par');
md=extrude(md,3,1.);
md=setflowequation(md,'SSA','all');
md.timestepping.time_step=0;
md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Thermal');

%Fields and tolerances to track changes
field_names     ={'Temperature','BasalforcingsGroundediceMeltingRate'};
field_tolerances={1e-13,6e-6};
field_values={...
	(md.results.ThermalSolution.Temperature),...
	(md.results.ThermalSolution.BasalforcingsGroundediceMeltingRate),...
	};
