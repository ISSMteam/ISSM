//Test Name: SquareShelfConstrainedStressHO
var md = new model();
triangle(md,square[0],180000.);
setmask(md,'all','');
parameterize(md);
md.extrude(md,3,2.);
setflowequation(md,'HO','all');
//md.cluster=generic('name',oshostname(),'np',3);
md.stressbalance.requested_outputs=['default','StressTensorxx','StressTensoryy','StressTensorzz','StressTensorxy','StressTensorxz','StressTensoryz'];
md=solve(md,'Stressbalance', 'checkconsistency', 'no');

//Fields and tolerances to track changes
field_names     =['Vx','Vy','Vz','Vel','Pressure',
	'StressTensorxx','StressTensoryy','StressTensorzz','StressTensorxy','StressTensorxz','StressTensoryz'];
field_tolerances=[1e-09,1e-09,1e-09,1e-09,1e-09,
	1e-09,1e-09,1e-09,1e-09,1e-09,1e-08];
field_values=[
	(md.results.StressbalanceSolution[0].Vx),
	(md.results.StressbalanceSolution[0].Vy),
	(md.results.StressbalanceSolution[0].Vz),
	(md.results.StressbalanceSolution[0].Vel),
	(md.results.StressbalanceSolution[0].Pressure),
	(md.results.StressbalanceSolution[0].StressTensorxx),
	(md.results.StressbalanceSolution[0].StressTensoryy),
	(md.results.StressbalanceSolution[0].StressTensorzz),
	(md.results.StressbalanceSolution[0].StressTensorxy),
	(md.results.StressbalanceSolution[0].StressTensorxz),
	(md.results.StressbalanceSolution[0].StressTensoryz),
	];
