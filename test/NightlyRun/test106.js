//Test Name: SquareShelfConstrainedMasstransp2dDG
var md = new model();
triangle(md,square[0],150000.);
meshconvert(md);
setmask(md,'all','');
parameterize(md);
setflowequation(md,'SSA','all');
//md.cluster=generic('name',oshostname(),'np',3);
md.masstransport.stabilization=3;
md.masstransport.spcthickness=md.geometry.thickness;
md=solve(md,'Masstransport');

//Fields and tolerances to track changes
field_names     =['Thickness'];
field_tolerances=[1e-13];
field_values=[
	(md.results.MasstransportSolution[0].Thickness),
	];
