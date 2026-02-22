%Test Name: PigShakti
md=triangle(model(),'../Exp/Pig.exp',1000.);
md=setmask(md,'../Exp/PigShelves.exp','../Exp/PigIslands.exp');
md=parameterize(md,'../Par/Pig.par');
md=setflowequation(md,'SSA','all');
md.timestepping.start_time = 0;
md.timestepping.time_step  = 1*3600/md.constants.yts;;
md.timestepping.final_time = 1/24/365;

% Turn on SHAKTI and turn off other transient processes for now
md.transient=deactivateall(md.transient);
md.transient.isstressbalance=0; % Solve for ice velocity 1, turn off 0
md.transient.ishydrology=1;

% HYDROLOGY SPECIFIC PARAMETERIZATION:
% Change hydrology class to SHAKTI model
md.hydrology=hydrologyshakti();

% Define distributed englacial input to the subglacial system (m/yr)
md.hydrology.englacial_input = 0.0*ones(md.mesh.numberofvertices,1);

% Define initial water head such that water pressure is 50% of ice overburden pressure
md.hydrology.head = 0.5*md.materials.rho_ice/md.materials.rho_freshwater*md.geometry.thickness + md.geometry.base;

% Initial subglacial gap height of 0.001m everywhere
md.hydrology.gap_height = 0.001*ones(md.mesh.numberofelements,1);

% Typical bed bump bump spacing
md.hydrology.bump_spacing = 1.0*ones(md.mesh.numberofelements,1);

% Typical bed bump height
md.hydrology.bump_height = 0.0*ones(md.mesh.numberofelements,1);

% Initial Reynolds number (start at Re=1000 everywhere)
md.hydrology.reynolds= 1000*ones(md.mesh.numberofelements,1);

% Deal with boundary conditions
md.hydrology.spchead = NaN(md.mesh.numberofvertices,1);
md.hydrology.spchead(md.mask.ocean_levelset<=0)=0;

md.hydrology.moulin_input = zeros(md.mesh.numberofvertices,1);
md.hydrology.neumannflux=zeros(md.mesh.numberofelements,1); 

% Friction
Neff = md.materials.rho_ice*md.constants.g*md.geometry.thickness-md.materials.rho_water*md.constants.g*(md.hydrology.head - md.geometry.base);
md.friction.effective_pressure=Neff;
md.friction.effective_pressure_limit(:)=0;
md.friction.coupling = 4;

md.cluster=generic('name',oshostname(),'np',8);
md=solve(md,'Transient');

%Fields and tolerances to track changes
field_names ={...
	'HydrologyHead1','HydrologyGapHeight1','EffectivePressure1'};
field_tolerances={...
	1e-13, 1e-13,...
	1e-13};
field_values={...
	md.results.TransientSolution(1).HydrologyHead, ...
	md.results.TransientSolution(1).HydrologyGapHeight,...
	md.results.TransientSolution(1).EffectivePressure};

