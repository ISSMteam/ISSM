//Test Name: SquareShelfTherTranForcTemp
var md = new model();
triangle(md,square[0],180000.);
setmask(md,'all','');
parameterize(md);
md.extrude(md,3,1.);
setflowequation(md,'SSA','all');
//md.cluster=generic('name',oshostname(),'np',3);
md.thermal.spctemperature=[[md.thermal.spctemperature, md.thermal.spctemperature+5., md.thermal.spctemperature+10., md.thermal.spctemperature+15.], [1.5, 2.5, 3.5, 4.]];
md.timestepping.time_step=1;
md.timestepping.final_time=4;
md.trans.isstressbalance=0;
md.trans.ismasstransport=0;
md.trans.issmb=1;
md.trans.isthermal=1;
md.trans.isgroundingline=0;
md=solve(md,'Transient');

//Fields and tolerances to track changes
field_names     =['Temperature1','BasalforcingsGroundediceMeltingRate1','Temperature2','BasalforcingsGroundediceMeltingRate2','Temperature3','BasalforcingsGroundediceMeltingRate3','Temperature4','BasalforcingsGroundediceMeltingRate4'];
field_tolerances=[1e-13,1e-6,1e-13,1e-6,1e-13,1e-6,1e-13,1e-6];
field_values=[
	(md.results.TransientSolution[0](1).Temperature),
	(md.results.TransientSolution[0](1).BasalforcingsGroundediceMeltingRate),
	(md.results.TransientSolution[0](2).Temperature),
	(md.results.TransientSolution[0](2).BasalforcingsGroundediceMeltingRate),
	(md.results.TransientSolution[0](3).Temperature),
	(md.results.TransientSolution[0](3).BasalforcingsGroundediceMeltingRate),
	(md.results.TransientSolution[0](4).Temperature),
	(md.results.TransientSolution[0](4).BasalforcingsGroundediceMeltingRate),
	];
