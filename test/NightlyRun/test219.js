//Test Name: SquareShelfStressSSAHOTiling
var md = new model();
triangle(md,square[0],150000.);
setmask(md,'all','');
parameterize(md);
md.extrude(md,3,2.);
setflowequation(md,'HO','../Exp/SquareHalfRight.exp','fill','SSA');
//md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Stressbalance');

//Fields and tolerances to track changes
field_names     =['Vx','Vy','Vz','Vel','Pressure'];
field_tolerances=[5e-09,5e-09,5e-09,5e-09,1e-13];
field_values=[
	(md.results.StressbalanceSolution[0].Vx),
	(md.results.StressbalanceSolution[0].Vy),
	(md.results.StressbalanceSolution[0].Vz),
	(md.results.StressbalanceSolution[0].Vel),
	(md.results.StressbalanceSolution[0].Pressure),
	];
