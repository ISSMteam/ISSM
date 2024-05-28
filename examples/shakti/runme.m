steps=[1:3];

if any(steps==1) 
	disp('	Step 1: Mesh');

	%Generate unstructured mesh on 1,000 m square with typical element edge length of 20 m
	md=triangle(model,'./outline.exp',20);

	save MoulinMesh md
end 

if any(steps==2) 
	disp('	Step 2: Parameterization');
	md=loadmodel('MoulinMesh');

	md=setmask(md,'','');

	% Run parameterization script to set up geometry, velocity, material properties, etc.
	md=parameterize(md,'moulin.par');

	% HYDROLOGY SPECIFIC PARAMETERIZATION:
	% Change hydrology class to Sommers' SHaKTI model
	md.hydrology=hydrologyshakti();

	% Define initial water head such that water pressure is 50% of ice overburden pressure
	md.hydrology.head = 0.5*md.materials.rho_ice/md.materials.rho_freshwater*md.geometry.thickness + md.geometry.base;

	% Initial subglacial gap height of 0.01m everywhere
	md.hydrology.gap_height = 0.01*ones(md.mesh.numberofelements,1);

	% Typical bed bump bump spacing (2m)
	md.hydrology.bump_spacing = 2*ones(md.mesh.numberofelements,1);

	% Typical bed bump height (0.1m)
	md.hydrology.bump_height = 0.1*ones(md.mesh.numberofelements,1);

	% Define distributed englacial input to the subglacial system (m/yr)
	% Change the value 0.0 to add distributed input
	md.hydrology.englacial_input = 0.0*ones(md.mesh.numberofvertices,1);

	% Initial Reynolds number (start at Re=1000 everywhere)
	md.hydrology.reynolds= 1000*ones(md.mesh.numberofelements,1);

	% Set up atmospheric pressure Type I boundary condition at left edge of
	% domain (outflow, i.e. h=zb at x=xmin)
	md.hydrology.spchead = NaN(md.mesh.numberofvertices,1);
	pos=find(md.mesh.vertexonboundary & md.mesh.x==min(md.mesh.x));
	md.hydrology.spchead(pos)=md.geometry.base(pos);

	save MoulinParam md;
end 

if any(steps==3) 
	disp('	Step 3: Solve!');
	md=loadmodel('MoulinParam');

	md.transient=deactivateall(md.transient);
	md.transient.ishydrology=1;

	% Specify that you want to run the model on your current computer
	% Change the number of processors according to your machine (here np=4)
	md.cluster=generic('np',2);

	% Define the time stepping scheme: run for 90 days with a time step of 1 hr
	md.timestepping.time_step=3600/md.constants.yts; % Time step (in years)
	md.timestepping.final_time=30/365;

	%Add one moulin with steady input at x=500, y=500
	[a,pos] = min(sqrt((md.mesh.x-500).^2+(md.mesh.y-500).^2));
	time=0:md.timestepping.time_step:md.timestepping.final_time;
	md.hydrology.moulin_input = zeros(md.mesh.numberofvertices+1,numel(time));
	md.hydrology.moulin_input(end,:)=time;
	md.hydrology.moulin_input(pos,:)=4;

	% Specify no-flux Type 2 boundary conditions on all edges (except
	% the Type 1 condition set at the outflow above)
	md.hydrology.neumannflux=zeros(md.mesh.numberofelements+1,numel(time));
	md.hydrology.neumannflux(end,:)=time;

	md.verbose.solution=1;
	md=solve(md,'Transient');

	save MoulinTransient md
end 
