//Test Name: SquareShelfTherTran
var md = new model();
triangle(md,square[0],180000.);
setmask(md,'all','');
parameterize(md);
md.extrude(md,3,1.);
setflowequation(md,'SSA','all');
//md.cluster=generic('name',oshostname(),'np',3);
md.transient.isstressbalance=0;
md.transient.ismasstransport=0;
md.transient.issmb=1;
md.transient.isthermal=1;
md.transient.isgroundingline=0;
md=solve(md,'Transient');

//Fields and tolerances to track changes
field_names     =['Temperature1','BasalforcingsGroundediceMeltingRate1','Temperature2','BasalforcingsGroundediceMeltingRate2','Temperature3','BasalforcingsGroundediceMeltingRate3'];
field_tolerances=[1e-13,1e-6,1e-13,1e-6,1e-13,1e-6];
field_values=[
	(md.results.TransientSolution[0](1).Temperature),
	(md.results.TransientSolution[0](1).BasalforcingsGroundediceMeltingRate),
	(md.results.TransientSolution[0](2).Temperature),
	(md.results.TransientSolution[0](2).BasalforcingsGroundediceMeltingRate),
	(md.results.TransientSolution[0](3).Temperature),
	(md.results.TransientSolution[0](3).BasalforcingsGroundediceMeltingRate),
	];
