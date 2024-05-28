//Test Name: SquareShelfStressSSA2d
var md = new model();
triangle(md,square[0],150000.);
setmask(md,'all','');
parameterize(md);
setflowequation(md,'SSA','all');
//md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Stressbalance');

//Fields and tolerances to track changes
field_names     =['Vx','Vy','Vel','Pressure'];
field_tolerances=[1e-13,1e-13,1e-13,1e-13];
field_values=[
	(md.results.StressbalanceSolution[0].Vx),
	(md.results.StressbalanceSolution[0].Vy),
	(md.results.StressbalanceSolution[0].Vel),
	(md.results.StressbalanceSolution[0].Pressure),
	];
