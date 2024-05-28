%Test Name: SquareSheetConstrainedTherTranNyeH2O
md=triangle(model(),'../Exp/Square.exp',180000.);
md=setmask(md,'','');
md=parameterize(md,'../Par/SquareSheetConstrained.par');
md=extrude(md,3,1.);
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',3);
md.verbose=verbose('convergence',true,'solution',true);
md.materials.rheology_law = 'NyeH2O';
md.materials.rheology_B=nye(md.initialization.temperature,2);

md.transient.isstressbalance=0;
md.transient.ismasstransport=0;
md.transient.issmb=1;
md.transient.isthermal=1;
md.transient.isgroundingline=0;
md=solve(md,'Transient');

%Fields and tolerances to track changes
field_names     ={'Temperature1','BasalforcingsGroundediceMeltingRate1',...
	'Temperature3','BasalforcingsGroundediceMeltingRate3'};
field_tolerances={1e-13,1e-13,1e-13,1e-13}; 
field_values={...
	(md.results.TransientSolution(1).Temperature),...
	(md.results.TransientSolution(1).BasalforcingsGroundediceMeltingRate),...
	(md.results.TransientSolution(3).Temperature),...
	(md.results.TransientSolution(3).BasalforcingsGroundediceMeltingRate),...
	};
