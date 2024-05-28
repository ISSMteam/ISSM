//Test Name: SquareShelfTranForcePos3d
var md = new model();
triangle(md,square[0],350000.);
setmask(md,'all','');
parameterize(md);
md.extrude(md,3,1.);
setflowequation(md,'SSA','all');
//md.cluster=generic('name',oshostname(),'np',3);

md.timestepping.time_step=1;
md.settings.output_frequency=1;
md.timestepping.final_time=4;

//Set up transient
smb=ones(md.mesh.numberofvertices,1);
for (var i = 0; i < smb.length; ++i) {
    smb[i][0] *= 3.6;
    smb[i].push(smb[i][0]*2);
}

md.smb.mass_balance=smb.slice();
md.smb.mass_balance[md.smb.mass_balance.length-1] = [1.5, 3.];
md.trans.isthermal=0;

md=solve(md,'Transient');

//Fields and tolerances to track changes
field_names     =['Vx1','Vy1','Vz1','Vel1','Pressure1','Bed1','Surface1','Thickness1','SmbMassBalance1', 
	'Vx2','Vy2','Vz2','Vel2','Pressure2','Bed2','Surface2','Thickness2','SmbMassBalance2', 
	'Vx3','Vy3','Vz3','Vel3','Pressure3','Bed3','Surface3','Thickness3','SmbMassBalance3', 
	'Vx4','Vy4','Vz4','Vel4','Pressure4','Bed4','Surface4','Thickness4','SmbMassbalance4'];
field_tolerances=[1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,
	1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,
	1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,
	1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13];
field_values=[
	(md.results.TransientSolution[0](1).Vx),
	(md.results.TransientSolution[0](1).Vy),
	(md.results.TransientSolution[0](1).Vz),
	(md.results.TransientSolution[0](1).Vel),
	(md.results.TransientSolution[0](1).Pressure),
	(md.results.TransientSolution[0](1).Base),
	(md.results.TransientSolution[0](1).Surface),
	(md.results.TransientSolution[0](1).Thickness),
	(md.results.TransientSolution[0](1).SmbMassBalance),
	(md.results.TransientSolution[0](2).Vx),
	(md.results.TransientSolution[0](2).Vy),
	(md.results.TransientSolution[0](2).Vz),
	(md.results.TransientSolution[0](2).Vel),
	(md.results.TransientSolution[0](2).Pressure),
	(md.results.TransientSolution[0](2).Base),
	(md.results.TransientSolution[0](2).Surface),
	(md.results.TransientSolution[0](2).Thickness),
	(md.results.TransientSolution[0](2).SmbMassBalance),
	(md.results.TransientSolution[0](3).Vx),
	(md.results.TransientSolution[0](3).Vy),
	(md.results.TransientSolution[0](3).Vz),
	(md.results.TransientSolution[0](3).Vel),
	(md.results.TransientSolution[0](3).Pressure),
	(md.results.TransientSolution[0](3).Base),
	(md.results.TransientSolution[0](3).Surface),
	(md.results.TransientSolution[0](3).Thickness),
	(md.results.TransientSolution[0](3).SmbMassBalance),
	(md.results.TransientSolution[0](4).Vx),
	(md.results.TransientSolution[0](4).Vy),
	(md.results.TransientSolution[0](4).Vz),
	(md.results.TransientSolution[0](4).Vel),
	(md.results.TransientSolution[0](4).Pressure),
	(md.results.TransientSolution[0](4).Base),
	(md.results.TransientSolution[0](4).Surface),
	(md.results.TransientSolution[0](4).Thickness),
	(md.results.TransientSolution[0](4).SmbMassBalance),
	];
