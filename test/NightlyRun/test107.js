//Test Name: SquareShelfConstrainedMasstransp3d
var md = new model();
triangle(md,square[0],150000.);
setmask(md,'all','');
parameterize(md);
setflowequation(md,'SSA','all');
md.extrude(md,5,3.);
//md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Masstransport');

//Fields and tolerances to track changes
field_names     =['Thickness'];
field_tolerances=[1e-13];
field_values=[
	(md.results.MasstransportSolution[0].Thickness),
	];
