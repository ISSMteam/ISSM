//Test Name: SquareShelfConstrainedTranSSA2d
var md = new model();
triangle(md,square[0],150000.);
setmask(md,'all','');
parameterize(md);
setflowequation(md,'SSA','all');
//md.cluster=generic('name',oshostname(),'np',3);
md.trans.requested_outputs=['IceVolume'];

md=solve(md,'transient');

//Fields and tolerances to track changes
field_names     =['Vx1','Vy1','Vel1','Pressure1','Bed1','Surface1','Thickness1','Volume1','Vx2','Vy2','Vel2','Pressure2','Bed2','Surface2','Thickness2','Volume2','Vx3','Vy3','Vel3','Pressure3','Bed3','Surface3','Thickness3','Volume3'];
field_tolerances=[1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,
						1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,
						1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13];
field_values=[
	(md.results.TransientSolution[0](1).Vx),
	(md.results.TransientSolution[0](1).Vy),
	(md.results.TransientSolution[0](1).Vel),
	(md.results.TransientSolution[0](1).Pressure),
	(md.results.TransientSolution[0](1).Base),
	(md.results.TransientSolution[0](1).Surface),
	(md.results.TransientSolution[0](1).Thickness),
	(md.results.TransientSolution[0](1).IceVolume),
	(md.results.TransientSolution[0](2).Vx),
	(md.results.TransientSolution[0](2).Vy),
	(md.results.TransientSolution[0](2).Vel),
	(md.results.TransientSolution[0](2).Pressure),
	(md.results.TransientSolution[0](2).Base),
	(md.results.TransientSolution[0](2).Surface),
	(md.results.TransientSolution[0](2).Thickness),
	(md.results.TransientSolution[0](2).IceVolume),
	(md.results.TransientSolution[0](3).Vx),
	(md.results.TransientSolution[0](3).Vy),
	(md.results.TransientSolution[0](3).Vel),
	(md.results.TransientSolution[0](3).Pressure),
	(md.results.TransientSolution[0](3).Base),
	(md.results.TransientSolution[0](3).Surface),
	(md.results.TransientSolution[0](3).Thickness),
	(md.results.TransientSolution[0](3).IceVolume),
	];
