%Test Name: SquareSheetHydrologyShakti
md=triangle(model(),'../Exp/Square.exp',50000.);
md.mesh.x = md.mesh.x/1000;
md.mesh.y = md.mesh.y/1000;
md=setmask(md,'','');
md=parameterize(md,'../Par/SquareSheetConstrained.par');
md.transient=deactivateall(md.transient);
md.transient.ishydrology=1;
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',2);

%Use hydrology coupled friction law
md.friction=frictionshakti(md.friction);

%Change hydrology class to Shakti' model
md.hydrology=hydrologyshakti();

%Change geometry
md.geometry.base = -.02*md.mesh.x + 20;
md.geometry.thickness = 300*ones(md.mesh.numberofvertices,1);
md.geometry.bed = md.geometry.base;
md.geometry.surface = md.geometry.bed+md.geometry.thickness;

%define the initial water head as being such that the water pressure is 50% of the ice overburden pressure
md.hydrology.head = 0.5*md.materials.rho_ice/md.materials.rho_freshwater*md.geometry.thickness + md.geometry.base;
md.hydrology.gap_height = 0.01*ones(md.mesh.numberofelements,1);
md.hydrology.bump_spacing = 2*ones(md.mesh.numberofelements,1);
md.hydrology.bump_height = 0.05*ones(md.mesh.numberofelements,1);
md.hydrology.englacial_input = 0.5*ones(md.mesh.numberofvertices,1);
md.hydrology.reynolds= 1000*ones(md.mesh.numberofelements,1);
md.hydrology.spchead = NaN(md.mesh.numberofvertices,1);
pos=find(md.mesh.vertexonboundary & md.mesh.x==1000);
md.hydrology.spchead(pos)=md.geometry.base(pos);

%Define velocity
md.initialization.vx = 10^-6*md.constants.yts*ones(md.mesh.numberofvertices,1);
md.initialization.vy = zeros(md.mesh.numberofvertices,1);

md.timestepping.time_step=3*3600/md.constants.yts;
md.timestepping.final_time=.5/365;
md.materials.rheology_B(:)= (5e-25)^(-1/3);

%Add one moulin and Neumann BC, varying in time
[a pos] = min(sqrt((md.mesh.x-500).^2+(md.mesh.y-500).^2));
time=0:md.timestepping.time_step:md.timestepping.final_time;
md.hydrology.moulin_input = zeros(md.mesh.numberofvertices+1,numel(time));
md.hydrology.moulin_input(end,:)=time;
md.hydrology.moulin_input(pos,:)=5*(1-sin(2*pi/(1/365)*time));
md.hydrology.neumannflux=zeros(md.mesh.numberofelements+1,numel(time));
md.hydrology.neumannflux(end,:)=time;
segx = md.mesh.x(md.mesh.segments(:,1)); segy=md.mesh.y(md.mesh.segments(:,1));
pos = md.mesh.segments(find(segx<1 & segy>400 & segy<600),3);
md.hydrology.neumannflux(pos,:)=repmat(0.05*(1-sin(2*pi/(1/365)*time)),numel(pos),1);

md=solve(md,'Transient');

%Fields and tolerances to track changes
field_names ={...
	'HydrologyHead1','HydrologyGapHeight1',...
	'HydrologyHead2','HydrologyGapHeight2',...
	'HydrologyHead3','HydrologyGapHeight3',...
	'HydrologyHead4','HydrologyGapHeight4'};
field_tolerances={...
	1e-13, 1e-13,...
	1e-13, 1e-13,...
	1e-13, 1e-13,...
	1e-13, 1e-12};
field_values={...
	md.results.TransientSolution(1).HydrologyHead, ...
	md.results.TransientSolution(1).HydrologyGapHeight,...
	md.results.TransientSolution(2).HydrologyHead, ...
	md.results.TransientSolution(2).HydrologyGapHeight,...
	md.results.TransientSolution(3).HydrologyHead, ...
	md.results.TransientSolution(3).HydrologyGapHeight,...
	md.results.TransientSolution(4).HydrologyHead, ...
	md.results.TransientSolution(4).HydrologyGapHeight};

