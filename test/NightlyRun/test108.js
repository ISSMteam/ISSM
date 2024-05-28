//Test Name: SquareShelfConstrainedTherStea
var md = new model();
triangle(md,square[0],180000.);
setmask(md,'all','');
parameterize(md);
md.extrude(md,3,1.);
setflowequation(md,'SSA','all');
md.timestepping.time_step=0;
//md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Thermal');

//Fields and tolerances to track changes
field_names     =['Temperature','BasalforcingsGroundediceMeltingRate'];
field_tolerances=[1e-13,1e-5];
field_values=[
	(md.results.ThermalSolution[0].Temperature),
	(md.results.ThermalSolution[0].BasalforcingsGroundediceMeltingRate),
	];
