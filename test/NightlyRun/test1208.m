%Test Name: EISMINTA
%EISMINT benchmark experiment A
numlayers=8;
resolution=50000.;

%To begin with the numerical model
md=triangle(model(),'../Exp/SquareEISMINT750000.exp',resolution);
md=setmask(md,'','');
md=parameterize(md,'../Par/RoundSheetEISMINT.par');

%We extrude the model to have a 3d model
md=extrude(md,numlayers,1.);
md=setflowequation(md,'SIA','all');

%Spc the nodes on the bed
pos=find(md.mesh.vertexonbase);
md.stressbalance.spcvx(pos)=0.;
md.stressbalance.spcvy(pos)=0.;
md.stressbalance.spcvz(pos)=0.;

%Adapt the time steps to the resolution
md.timestepping.time_step=15.;
md.settings.output_frequency=500;
md.timestepping.final_time=30000.;
md.masstransport.stabilization=1;
md.thermal.stabilization=1;

%Now we can solve the problem 
md.cluster=generic('name',oshostname(),'np',8);
md=solve(md,'Transient');

%Fields and tolerances to track changes
field_names     ={'Vx','Vy','Vz','Vel','Pressure','Thickness','Base','Surface','Temperature','BasalforcingsGroundediceMeltingRate'};
field_tolerances={1e-08,1e-08,1e-07,1e-08,1e-08,1e-08,1e-08,1e-08,1e-07,3e-07};
field_values={...
	(md.results.TransientSolution(end).Vx),...
	(md.results.TransientSolution(end).Vy),...
	(md.results.TransientSolution(end).Vz),...
	(md.results.TransientSolution(end).Vel),...
	(md.results.TransientSolution(end).Pressure),...
	(md.results.TransientSolution(end).Thickness),...
	(md.results.TransientSolution(end).Base),...
	(md.results.TransientSolution(end).Surface),...
	(md.results.TransientSolution(end).Temperature),...
	(md.results.TransientSolution(end).BasalforcingsGroundediceMeltingRate),...
	};
