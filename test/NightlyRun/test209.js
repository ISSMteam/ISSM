//Test Name: SquareShelfTranSSA3d
var md = new model();
triangle(md,square[0],180000.);
setmask(md,'all','');
parameterize(md);
md.extrude(md,3,1.);
setflowequation(md,'SSA','all');
//md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Transient');

//Fields and tolerances to track changes
field_names     =['Vx1','Vy1','Vz1','Vel1','Pressure1','Bed1','Surface1','Thickness1','Temperature1','BasalforcingsGroundediceMeltingRate1', 
	'Vx2','Vy2','Vz2','Vel2','Pressure2','Bed2','Surface2','Thickness2','Temperature2','BasalforcingsGroundediceMeltingRate2', 
	'Vx3','Vy3','Vz3','Vel3','Pressure3','Bed3','Surface3','Thickness3','Temperature3','BasalforcingsGroundediceMeltingRate3'];
field_tolerances=[
	1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-8,
	1e-13,1e-13,3e-9,1e-13,3.8e-11,3.85e-11,3.85e-11,3.85e-11,1e-13,2e-8,
	5e-11,6e-12,1e-08,6e-12,3.9e-11,1e-10,1e-10,1e-10,9e-12,5e-8];
field_values=[
	(md.results.TransientSolution[0](1).Vx),
	(md.results.TransientSolution[0](1).Vy),
	(md.results.TransientSolution[0](1).Vz),
	(md.results.TransientSolution[0](1).Vel),
	(md.results.TransientSolution[0](1).Pressure),
	(md.results.TransientSolution[0](1).Base),
	(md.results.TransientSolution[0](1).Surface),
	(md.results.TransientSolution[0](1).Thickness),
	(md.results.TransientSolution[0](1).Temperature),
	(md.results.TransientSolution[0](1).BasalforcingsGroundediceMeltingRate),
	(md.results.TransientSolution[0](2).Vx),
	(md.results.TransientSolution[0](2).Vy),
	(md.results.TransientSolution[0](2).Vz),
	(md.results.TransientSolution[0](2).Vel),
	(md.results.TransientSolution[0](2).Pressure),
	(md.results.TransientSolution[0](2).Base),
	(md.results.TransientSolution[0](2).Surface),
	(md.results.TransientSolution[0](2).Thickness),
	(md.results.TransientSolution[0](2).Temperature),
	(md.results.TransientSolution[0](2).BasalforcingsGroundediceMeltingRate),
	(md.results.TransientSolution[0](3).Vx),
	(md.results.TransientSolution[0](3).Vy),
	(md.results.TransientSolution[0](3).Vz),
	(md.results.TransientSolution[0](3).Vel),
	(md.results.TransientSolution[0](3).Pressure),
	(md.results.TransientSolution[0](3).Base),
	(md.results.TransientSolution[0](3).Surface),
	(md.results.TransientSolution[0](3).Thickness),
	(md.results.TransientSolution[0](3).Temperature),
	(md.results.TransientSolution[0](3).BasalforcingsGroundediceMeltingRate),
	];
