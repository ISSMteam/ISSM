%Test Name: SquareSheetShelfTherStea
md=triangle(model(),'../Exp/Square.exp',150000.);
md=setmask(md,'../Exp/SquareShelf.exp','');
md=parameterize(md,'../Par/SquareSheetShelf.par');
md=extrude(md,4,1.);
md=setflowequation(md,'HO','all');
md.cluster=generic('name',oshostname(),'np',3);
md.timestepping.time_step=0.;
md=solve(md,'Thermal');

%Fields and tolerances to track changes
field_names     ={'Temperature','BasalforcingsGroundediceMeltingRate'};
field_tolerances={1e-13,1e-5};
field_values={...
	(md.results.ThermalSolution.Temperature),...
	(md.results.ThermalSolution.BasalforcingsGroundediceMeltingRate),...
	};
