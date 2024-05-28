%Test Name: PigTherTranSUPG
md=triangle(model(),'../Exp/Pig.exp',30000.);
md=setmask(md,'../Exp/PigShelves.exp','../Exp/PigIslands.exp');
md=parameterize(md,'../Par/Pig.par');
md=extrude(md,3,1.);
md=setflowequation(md,'HO','all');
md.thermal.stabilization=2;
md.cluster=generic('name',oshostname(),'np',3);
md.transient.isstressbalance=0;
md.transient.ismasstransport=0;
md.transient.issmb=1;
md.transient.isthermal=1;
md.transient.isgroundingline=0;
md=solve(md,'Transient');

%Fields and tolerances to track changes
field_names     ={'Temperature1','BasalforcingsGroundediceMeltingRate1', ...
				      'Temperature2','BasalforcingsGroundediceMeltingRate2'};
field_tolerances={1e-13,1e-8,1e-13,5e-8};
field_values={...
	(md.results.TransientSolution(1).Temperature),...
	(md.results.TransientSolution(1).BasalforcingsGroundediceMeltingRate),...
	(md.results.TransientSolution(2).Temperature),...
	(md.results.TransientSolution(2).BasalforcingsGroundediceMeltingRate),...
	};
