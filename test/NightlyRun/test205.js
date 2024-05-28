//Test Name: SquareShelfStressMHOPenalties
var md = new model();
triangle(md,square[0],150000.);
setmask(md,'all','');
parameterize(md);
md.extrude(md,3,2.);
setflowequation(md,'HO','../Exp/SquareHalfRight.exp','fill','SSA','coupling','penalties');
md.settings.solver_residue_threshold = 1.e-4;
//md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Stressbalance');

//Fields and tolerances to track changes
field_names     =['Vx','Vy','Vz','Vel','Pressure'];
field_tolerances=[2e-05,2e-05,1e-05,1e-05,1e-05];
field_values=[
	(md.results.StressbalanceSolution[0].Vx),
	(md.results.StressbalanceSolution[0].Vy),
	(md.results.StressbalanceSolution[0].Vz),
	(md.results.StressbalanceSolution[0].Vel),
	(md.results.StressbalanceSolution[0].Pressure),
	];
