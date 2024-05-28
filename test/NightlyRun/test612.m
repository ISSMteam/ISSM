%Test Name: 79NorthPISMhydrofriction
md=triangle(model(),'../Exp/79North.exp',10000.);
md=setmask(md,'../Exp/79NorthShelf.exp','');
md=parameterize(md,'../Par/79North.par');
md=setflowequation(md,'SSA','all');

%Hydrology
md.hydrology = hydrologypism();
md.hydrology.drainage_rate  = 0.001*ones(md.mesh.numberofvertices,1);
md.hydrology.watercolumn_max=2*ones(md.mesh.numberofvertices,1);
md.initialization.pressure=zeros(md.mesh.numberofvertices,1);
md.initialization.watercolumn=zeros(md.mesh.numberofvertices,1);
md.basalforcings.groundedice_melting_rate = [1:md.mesh.numberofvertices]';

%Friction
md.friction=frictionpism();
md.friction.till_friction_angle=30*ones(md.mesh.numberofvertices,1);
md.friction.sediment_compressibility_coefficient=0.12*ones(md.mesh.numberofvertices,1);

md.transient.ishydrology        = 1;
md.transient.issmb              = 0;
md.transient.ismasstransport    = 0;
md.transient.isstressbalance    = 1;
md.transient.isthermal          = 0;
md.transient.isgroundingline    = 0;

md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Transient');

%Plot to check result
% plotmodel(md,'data',md.results.TransientSolution(3).Vx - md.results.TransientSolution(2).Vx,'caxis',[-1 1])
% plotmodel(md,'data',md.results.TransientSolution(3).Vy - md.results.TransientSolution(2).Vy,'caxis',[-1 1])
% plotmodel(md,'data',md.results.TransientSolution(3).Vel - md.results.TransientSolution(2).Vel);

%Fields and tolerances to track changes
field_names     ={'Vx1','Vx2','Vx3','Vy1','Vy2','Vy3','Vel1','Vel2','Vel3'};
field_tolerances={2e-11,2e-11,2e-11,2e-11,2e-11,2e-11,2e-11,2e-11,2e-11};
field_values={...
	(md.results.TransientSolution(1).Vx),...
	(md.results.TransientSolution(2).Vx),...
	(md.results.TransientSolution(3).Vx),...
	(md.results.TransientSolution(1).Vy),...
	(md.results.TransientSolution(2).Vy),...
	(md.results.TransientSolution(3).Vy),...
	(md.results.TransientSolution(1).Vel),...
	(md.results.TransientSolution(2).Vel),...
	(md.results.TransientSolution(3).Vel),...
	};
