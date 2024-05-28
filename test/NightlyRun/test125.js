//Test Name: SquareShelfConstrainedRestartTranSSA2d
var md = new model();
triangle(md,square[0],150000.);
setmask(md,'all','');
parameterize(md);
setflowequation(md,'SSA','all');
//md.cluster=generic('name',oshostname(),'np',1);
md.transient.requested_outputs=['IceVolume','TotalSmb'];

//md.verbose=verbose('solution',true);
md.settings.checkpoint_frequency=5;

// time steps and resolution
md.timestepping.final_time=8;

md=solve(md,'TransientSolution');
md2=solve(md,'TransientSolution','restart',1);

//Fields and tolerances to track changes
field_names     =['Vx1','Vy1','Vel1','TotalSmb1','Bed1','Surface1','Thickness1','Volume1','Vx2','Vy2','Vel2','TotalSmb2','Bed2','Surface2','Thickness2','Volume2','Vx3','Vy3','Vel3','TotalSmb3','Bed3','Surface3','Thickness3','Volume3'];
field_tolerances=[1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,
						1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,
						1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13];
field_values=[
	(md.results.TransientSolution(6).Vx)-(md2.results.TransientSolution[0](1).Vx),
	(md.results.TransientSolution(6).Vy)-(md2.results.TransientSolution[0](1).Vy),
	(md.results.TransientSolution(6).Vel)-(md2.results.TransientSolution[0](1).Vel),
	(md.results.TransientSolution(6).TotalSmb)-(md2.results.TransientSolution[0](1).TotalSmb),
	(md.results.TransientSolution(6).Base)-(md2.results.TransientSolution[0](1).Base),
	(md.results.TransientSolution(6).Surface)-(md2.results.TransientSolution[0](1).Surface),
	(md.results.TransientSolution(6).Thickness)-(md2.results.TransientSolution[0](1).Thickness),
	(md.results.TransientSolution(6).IceVolume)-(md2.results.TransientSolution[0](1).IceVolume),
	(md.results.TransientSolution(7).Vx)-(md2.results.TransientSolution[0](2).Vx),
	(md.results.TransientSolution(7).Vy)-(md2.results.TransientSolution[0](2).Vy),
	(md.results.TransientSolution(7).Vel)-(md2.results.TransientSolution[0](2).Vel),
	(md.results.TransientSolution(7).TotalSmb)-(md2.results.TransientSolution[0](2).TotalSmb),
	(md.results.TransientSolution(7).Base)-(md2.results.TransientSolution[0](2).Base),
	(md.results.TransientSolution(7).Surface)-(md2.results.TransientSolution[0](2).Surface),
	(md.results.TransientSolution(7).Thickness)-(md2.results.TransientSolution[0](2).Thickness),
	(md.results.TransientSolution(7).IceVolume)-(md2.results.TransientSolution[0](2).IceVolume),
	(md.results.TransientSolution(8).Vx)-(md2.results.TransientSolution[0](3).Vx),
	(md.results.TransientSolution(8).Vy)-(md2.results.TransientSolution[0](3).Vy),
	(md.results.TransientSolution(8).Vel)-(md2.results.TransientSolution[0](3).Vel),
	(md.results.TransientSolution(8).TotalSmb)-(md2.results.TransientSolution[0](3).TotalSmb),
	(md.results.TransientSolution(8).Base)-(md2.results.TransientSolution[0](3).Base),
	(md.results.TransientSolution(8).Surface)-(md2.results.TransientSolution[0](3).Surface),
	(md.results.TransientSolution(8).Thickness)-(md2.results.TransientSolution[0](3).Thickness),
	(md.results.TransientSolution(8).IceVolume)-(md2.results.TransientSolution[0](3).IceVolume),
	];
