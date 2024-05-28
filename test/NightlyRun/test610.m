%Test Name: 79NorthPISMhydro
md=triangle(model(),'../Exp/79North.exp',10000.);
md=setmask(md,'../Exp/79NorthShelf.exp','');
md=parameterize(md,'../Par/79North.par');
md=setflowequation(md,'SSA','all');

%Hydrology
md.hydrology = hydrologypism();
md.hydrology.drainage_rate      = 0.001*ones(md.mesh.numberofvertices,1);
md.hydrology.watercolumn_max=20*ones(md.mesh.numberofvertices,1);
md.initialization.pressure=zeros(md.mesh.numberofvertices,1);
md.initialization.watercolumn=zeros(md.mesh.numberofvertices,1);
md.basalforcings.groundedice_melting_rate = [1:md.mesh.numberofvertices]';
md.transient.ishydrology        = 1;
md.transient.issmb              = 0;
md.transient.ismasstransport    = 0;
md.transient.isstressbalance    = 0;
md.transient.isthermal          = 0;
md.transient.isgroundingline    = 0;

md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Transient');

%Plot to check result
% plotmodel(md,'data',md.results.TransientSolution(3).Watercolumn - 3*(md.materials.rho_freshwater/md.materials.rho_ice*[1:md.mesh.numberofvertices]' -1))

%Fields and tolerances to track changes
field_names     ={'WaterColumn1','WaterColumn2','WaterColumn3'};
field_tolerances={1e-12,1e-12,1e-12};
field_values={...
	(md.results.TransientSolution(1).Watercolumn),...
	(md.results.TransientSolution(2).Watercolumn),...
	(md.results.TransientSolution(3).Watercolumn),...
	};
