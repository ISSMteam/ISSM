%Test Name: SquareShelfTherTranForcTemp
md=triangle(model(),'../Exp/Square.exp',180000.);
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareShelf.par');
md=extrude(md,3,1.);
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',3);
md.thermal.spctemperature=[md.thermal.spctemperature, md.thermal.spctemperature+5., md.thermal.spctemperature+10., md.thermal.spctemperature+15.; 1.5 2.5 3.5 4.];
md.timestepping.time_step=1;
md.timestepping.final_time=4;
md.transient.isstressbalance=0;
md.transient.ismasstransport=0;
md.transient.issmb=1;
md.transient.isthermal=1;
md.transient.isgroundingline=0;
md=solve(md,'Transient');

%Fields and tolerances to track changes
field_names     ={'Temperature1','BasalforcingsGroundediceMeltingRate1','Temperature2','BasalforcingsGroundediceMeltingRate2','Temperature3','BasalforcingsGroundediceMeltingRate3','Temperature4','BasalforcingsGroundediceMeltingRate4'};
field_tolerances={1e-13,1e-6,1e-13,1e-6,1e-13,1e-6,1e-13,1e-6};
field_values={...
	(md.results.TransientSolution(1).Temperature),...
	(md.results.TransientSolution(1).BasalforcingsGroundediceMeltingRate),...
	(md.results.TransientSolution(2).Temperature),...
	(md.results.TransientSolution(2).BasalforcingsGroundediceMeltingRate),...
	(md.results.TransientSolution(3).Temperature),...
	(md.results.TransientSolution(3).BasalforcingsGroundediceMeltingRate),...
	(md.results.TransientSolution(4).Temperature),...
	(md.results.TransientSolution(4).BasalforcingsGroundediceMeltingRate),...
	};
