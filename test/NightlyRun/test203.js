//Test Name: SquareShelfStressHO
var md = new model();
triangle(md,square[0],180000.);
setmask(md,'all','');
parameterize(md);
md.extrude(md,3,2.);
setflowequation(md,'HO','all');
//md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Stressbalance');

//Fields and tolerances to track changes
field_names     =['Vx','Vy','Vz','Vel','Pressure'];
field_tolerances=[2e-09,2e-09,1e-09,2e-09,1e-09];
field_values=[
	(md.results.StressbalanceSolution[0].Vx),
	(md.results.StressbalanceSolution[0].Vy),
	(md.results.StressbalanceSolution[0].Vz),
	(md.results.StressbalanceSolution[0].Vel),
	(md.results.StressbalanceSolution[0].Pressure),
	];
