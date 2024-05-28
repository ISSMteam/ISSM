//Test Name: SquareShelfStressSSA3d
var md = new model();
triangle(md,square[0],180000.);
setmask(md,'all','');
parameterize(md);
md.extrude(md,3,2.);
setflowequation(md,'SSA','all');
//md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Stressbalance');

//Fields and tolerances to track changes
field_names     =['Vx','Vy','Vz','Vel','Pressure'];
field_tolerances=[1e-13,1e-13,1e-13,1e-13,1e-13];
field_values=[
	(md.results.StressbalanceSolution[0].Vx),
	(md.results.StressbalanceSolution[0].Vy),
	(md.results.StressbalanceSolution[0].Vz),
	(md.results.StressbalanceSolution[0].Vel),
	(md.results.StressbalanceSolution[0].Pressure),
	];
