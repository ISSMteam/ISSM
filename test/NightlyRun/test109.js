//Test Name: SquareShelfConstrainedTherTran
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
field_names     =['Temperature','BasalforcingsGroundediceMeltingRate'];
field_tolerances=[1e-13,1e-13];
field_values=[
	(md.results.TransientSolution[0](1).Temperature),
	(md.results.TransientSolution[0](1).BasalforcingsGroundediceMeltingRate),
	];
