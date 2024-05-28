//Test Name: SquareShelfStressSSAFS3dTiling
var md = new model();
triangle(md,square[0],120000.);
setmask(md,'all','');
parameterize(md);
md.extrude(md,2,1.);
setflowequation(md,'FS','../Exp/SquareHalfRight.exp','fill','SSA');
//md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Stressbalance');

//Fields and tolerances to track changes
field_names     =['Vx','Vy','Vz','Vel'];
field_tolerances=[1e-09,1e-09,5e-06,1e-09];
field_values=[
	(md.results.StressbalanceSolution[0].Vx),
	(md.results.StressbalanceSolution[0].Vy),
	(md.results.StressbalanceSolution[0].Vz),
	(md.results.StressbalanceSolution[0].Vel)
];
