//Test Name: SquareShelfStressFS
var md = new model();
triangle(md,square[0],180000.);
setmask(md,'all','');
parameterize(md);
md.extrude(md,3,2.);
setflowequation(md,'FS','all');
//md.cluster=generic('name',oshostname(),'np',1);
md.stressbalance.shelf_dampening=1;
md.timestepping.time_step=0;
md1=solve(md,'Stressbalance');
md.stressbalance.shelf_dampening=0;
md=solve(md,'Stressbalance');

//Fields and tolerances to track changes
field_names     =['Vx','Vy','Vz','Vel','Pressure','Vx_damp','Vy_damp','Vz_damp','Vel_damp','Pressure_damp'];
field_tolerances=[1e-08,1e-08,8e-06,1e-08,1e-08,1e-08,1e-08,2e-07,1e-08,1e-08];
field_values=[
	(md.results.StressbalanceSolution[0].Vx),
	(md.results.StressbalanceSolution[0].Vy),
	(md.results.StressbalanceSolution[0].Vz),
	(md.results.StressbalanceSolution[0].Vel),
	(md.results.StressbalanceSolution[0].Pressure),
	(md1.results.StressbalanceSolution.Vx),
	(md1.results.StressbalanceSolution.Vy),
	(md1.results.StressbalanceSolution.Vz),
	(md1.results.StressbalanceSolution.Vel),
	(md1.results.StressbalanceSolution.Pressure),
	];
