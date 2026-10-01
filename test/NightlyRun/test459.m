%Test Name: SquareSheetShelfHydrologyGlaDSTransitionModel

%create model:
md=triangle(model(),'../Exp/Square.exp',50000.);
md=setmask(md,'../Exp/SquareShelf.exp','');
md=parameterize(md,'../Par/SquareSheetShelf.par');
md.mesh.x = md.mesh.x/100;
md.mesh.y = md.mesh.y/100;
md.miscellaneous.name='testChannels';

%miscellaneous
md=setflowequation(md,'SSA','all');

%Some constants
md.constants.g=9.8;
md.materials.rho_ice=910;

%Define initial conditions
md.initialization.vx = 1.0e-6*md.constants.yts*ones(md.mesh.numberofvertices,1);
md.initialization.vy = zeros(md.mesh.numberofvertices,1);
md.initialization.temperature=(273.-20.)*ones(md.mesh.numberofvertices,1);
md.initialization.watercolumn=0.03*ones(md.mesh.numberofvertices,1);
md.initialization.hydraulic_potential = md.materials.rho_ice*md.constants.g*md.geometry.thickness;

%Materials
md.materials.rheology_B = (5e-25)^(-1/3) * ones(md.mesh.numberofvertices,1);
md.materials.rheology_n = 3.*ones(md.mesh.numberofelements,1);

%Friction
md.friction.coefficient=zeros(md.mesh.numberofvertices,1);
md.friction.p=ones(md.mesh.numberofelements,1);
md.friction.q=ones(md.mesh.numberofelements,1);

%Boundary conditions:
md=SetIceSheetBC(md);

md.transient=deactivateall(md.transient);
md.transient.ishydrology=1;
md.transient.isgroundingline = 1;

% Set numerical conditions
md.timestepping.time_step=.1/365;
md.timestepping.final_time=.4/365;

%Change hydrology class to Glads model
md.hydrology=hydrologyglads();
md.hydrology.maxiter=10; %Make sure it runs quickly...
md.hydrology.ischannels=1;
md.hydrology.istransition = 1;
md.hydrology.omega = 1/2000.;
md.hydrology.sheet_alpha = 3./2.;
md.hydrology.sheet_beta = 3./2.;
md.hydrology.englacial_void_ratio=1e-5;
md.hydrology.moulin_input=zeros(md.mesh.numberofvertices,1);
md.hydrology.neumannflux=zeros(md.mesh.numberofelements,1);
md.hydrology.bump_height = 1e-1 * ones(md.mesh.numberofvertices,1);
md.hydrology.sheet_conductivity= 1e-2 * ones(md.mesh.numberofvertices,1);
md.hydrology.channel_conductivity= 5.e-2 * ones(md.mesh.numberofvertices,1);
md.hydrology.rheology_B_base = cuffey(273.15 - 2)*ones(md.mesh.numberofvertices,1);

% BCs for hydrology
pos=find(md.mesh.y<=100 & md.mesh.vertexonboundary);
md.hydrology.spcphi=NaN(md.mesh.numberofvertices,1);
md.hydrology.spcphi(pos) = md.materials.rho_ice * md.constants.g * md.geometry.thickness(pos);

md.cluster=generic('np',2);
md.hydrology.requested_outputs = {'default', 'TotalHydrologyGroundinglineDischarge'};
md=solve(md,'Transient'); %or 'tr'

%Fields and tolerances to track changes
field_names ={...
	'HydrologySheetThickness1','HydraulicPotential1','ChannelArea1','TotalHydrologyGroundinglineDischarge'...
	'HydrologySheetThickness2','HydraulicPotential2','ChannelArea2','TotalHydrologyGroundinglineDischarge'...
	'HydrologySheetThickness3','HydraulicPotential3','ChannelArea3','TotalHydrologyGroundinglineDischarge'...
	'HydrologySheetThickness4','HydraulicPotential4','ChannelArea4','TotalHydrologyGroundinglineDischarge'};
field_tolerances={...
	9e-14, 3e-13, 3e-12, 1e-13,...
	2e-13, 3e-13, 3e-12, 1e-13,...
	3e-13, 3e-13, 3e-12, 1e-13,...
	3e-13, 2e-13, 3e-12, 1e-13};
field_values={...
	md.results.TransientSolution(1).HydrologySheetThickness, ...
	md.results.TransientSolution(1).HydraulicPotential,...
	md.results.TransientSolution(1).ChannelArea,...
	md.results.TransientSolution(1).TotalHydrologyGroundinglineDischarge,...
	md.results.TransientSolution(2).HydrologySheetThickness, ...
	md.results.TransientSolution(2).HydraulicPotential,...
	md.results.TransientSolution(2).ChannelArea,...
	md.results.TransientSolution(2).TotalHydrologyGroundinglineDischarge,...
	md.results.TransientSolution(3).HydrologySheetThickness, ...
	md.results.TransientSolution(3).HydraulicPotential,...
	md.results.TransientSolution(3).ChannelArea,...
	md.results.TransientSolution(3).TotalHydrologyGroundinglineDischarge,...
	md.results.TransientSolution(4).HydrologySheetThickness, ...
	md.results.TransientSolution(4).HydraulicPotential,...
	md.results.TransientSolution(4).ChannelArea,...
	md.results.TransientSolution(4).TotalHydrologyGroundinglineDischarge};

return

%plot channels and grounding line (MATLAB ONLY)
plotchannels(md,md.results.TransientSolution(end).ChannelArea);
contours=isoline(md, md.results.TransientSolution(end).MaskOceanLevelset,'output','matrix');
plot(contours(:,1),contours(:,2),'-r');
