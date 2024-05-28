//Test Name: SquareShelfStressSSA2dRift
var md = new model();
triangle(md,rifts[0],50000.);
meshprocessrifts(md,'../Exp/Square.exp');
setmask(md,'all','');
parameterize(md);
setflowequation(md,'SSA','all');
//md.cluster=generic('name',oshostname(),'np',3);

//rift settings
md.rifts.riftstruct.fill='Melange';
md.rifts.riftstruct.fraction=0;
md.stressbalance.rift_penalty_lock=2;
md.stressbalance.rift_penalty_threshold=0;
md.rifts.riftstruct.fractionincrement=.1;
md=solve(md,'Stressbalance');

//Fields and tolerances to track changes
field_names     =['Vx','Vy','Vel','Pressure'];
field_tolerances=[9e-7,7e-8,9e-8,2e-11];
field_values=[
	(md.results.StressbalanceSolution[0].Vx),
	(md.results.StressbalanceSolution[0].Vy),
	(md.results.StressbalanceSolution[0].Vel),
	(md.results.StressbalanceSolution[0].Pressure),
	];
