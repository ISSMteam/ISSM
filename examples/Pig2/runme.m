steps=[1];

if any(steps==1) %Mesh Generation #1 

	%Mesh parameters
	domain =['./DomainOutline.exp'];
	hinit=5000;		% element size for the initial mesh
	hmax=40000;		% maximum element size of the final mesh
	hmin=4000;		% minimum element size of the final mesh
	gradation=1.7;	% maximum size ratio between two neighboring elements
	err=8;			% maximum error between interpolated and control field

	% Generate an initial uniform mesh (resolution = hinit m)
	md=bamg(model,'domain',domain,'hmax',hinit);

	% Get necessary data to build up the velocity grid
	nsidc_vel	='../Data/Antarctica_ice_velocity.nc';
	xmin		= strsplit(ncreadatt(nsidc_vel,'/','xmin'));	xmin	= str2num(xmin{2});
	ymax		= strsplit(ncreadatt(nsidc_vel,'/','ymax'));	ymax	= str2num(ymax{2});
	spacing		= strsplit(ncreadatt(nsidc_vel,'/','spacing'));	spacing	= str2num(spacing{2});
	nx			= double(ncreadatt(nsidc_vel,'/','nx'));
	ny			= double(ncreadatt(nsidc_vel,'/','ny'));
	vx			= double(ncread(nsidc_vel,'vx'));
	vy			= double(ncread(nsidc_vel,'vy'));

	% Build the coordinates
	x=xmin+(0:1:nx)'*spacing;
	y=(ymax-ny*spacing)+(0:1:ny)'*spacing;

	% Interpolate velocities onto coarse mesh
	vx_obs=InterpFromGridToMesh(x,y,flipud(vx'),md.mesh.x,md.mesh.y,0);
	vy_obs=InterpFromGridToMesh(x,y,flipud(vy'),md.mesh.x,md.mesh.y,0);
	vel_obs=sqrt(vx_obs.^2+vy_obs.^2);
	clear vx vy x y;

	% Adapt the mesh to minimize error in velocity interpolation
	md=bamg(md,'hmax',hmax,'hmin',hmin,'gradation',gradation,'field',vel_obs,'err',err);

	% Plot and save model
	plotmodel(md,'data','mesh')
	save ./Models/PIG_Mesh_generation md;
end 

if any(steps==2) %Masks #2 
	md = loadmodel('./Models/PIG_Mesh_generation');

	% Load SeaRISe dataset for Antarctica  http://websrv.cs.umt.edu/isis/index.php/Present_Day_Antarctica
	searise='../Data/Antarctica_5km_withshelves_v0.75.nc';

	%read thickness mask from SeaRISE
	x1=double(ncread(searise,'x1'));
	y1=double(ncread(searise,'y1'));
	thkmask=double(ncread(searise,'thkmask'));

	%interpolate onto our mesh vertices
	groundedice=double(InterpFromGridToMesh(x1,y1,thkmask',md.mesh.x,md.mesh.y,0));
	groundedice(groundedice<=0)=-1;
	clear thkmask;

	%fill in the md.mask structure
	md.mask.ocean_levelset=groundedice; %ice is grounded for mask equal one
	md.mask.ice_levelset=-1*ones(md.mesh.numberofvertices,1);%ice is present when negatvie

	plotmodel(md,'data',md.mask.ocean_levelset,'title','grounded/floating','data',md.mask.ice_levelset,'title','ice/no-ice')

	save ./Models/PIG_SetMask md;
end 

if any(steps==3) %Parameterization #3 
	md = loadmodel('./Models/PIG_SetMask');
	md = setflowequation(md,'SSA','all');
	md = parameterize(md,'./Pig.par');

	save ./Models/PIG_Parameterization md;
end 

if any(steps==4) %Rheology B inversion 
	md = loadmodel('./Models/PIG_Parameterization');

	% Control general
	md.inversion.iscontrol=1;
	md.inversion.maxsteps=40;
	md.inversion.maxiter=40;
	md.inversion.dxmin=0.1;
	md.inversion.gttol=1.0e-6;
	md.verbose=verbose('control',true);

	% Cost functions
	md.inversion.cost_functions=[101 103 502];
	md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,3);
	md.inversion.cost_functions_coefficients(:,1)=1000;
	md.inversion.cost_functions_coefficients(:,2)=1;
	md.inversion.cost_functions_coefficients(:,3)=1.e-16;

	% Controls
	md.inversion.control_parameters={'MaterialsRheologyBbar'};
	md.inversion.min_parameters=md.materials.rheology_B;
	md.inversion.max_parameters=md.materials.rheology_B;
	pos = find(md.mask.ocean_levelset<0);
	md.inversion.min_parameters(pos) = cuffey(273);
	md.inversion.max_parameters(pos) = cuffey(200);

	% Additional parameters
	md.stressbalance.restol=0.01;
	md.stressbalance.reltol=0.1;
	md.stressbalance.abstol=NaN;

	% Solve
	md.cluster=generic('name',oshostname,'np',2);
	mds=extract(md,md.mask.ocean_levelset<0);
	mds=solve(mds,'Stressbalance');

	% Update model rheology_B accordingly
	md.materials.rheology_B(mds.mesh.extractedvertices)=mds.results.StressbalanceSolution.MaterialsRheologyBbar;
	plotmodel(md,'data',md.materials.rheology_B)

	% Save model
	save ./Models/PIG_Control_B md;
end 

if any(steps==5) %drag inversion 
	md = loadmodel('./Models/PIG_Control_B');

	% Cost functions
	md.inversion.cost_functions=[101 103 501];
	md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,3);
	md.inversion.cost_functions_coefficients(:,1)=2000;
	md.inversion.cost_functions_coefficients(:,2)=1;
	md.inversion.cost_functions_coefficients(:,3)=8e-7;

	% Controls
	md.inversion.control_parameters={'FrictionCoefficient'};
	md.inversion.min_parameters=1*ones(md.mesh.numberofvertices,1);
	md.inversion.max_parameters=200*ones(md.mesh.numberofvertices,1);

	% Solve
	md=solve(md,'Stressbalance');

	% Update model friction fields accordingly
	md.friction.coefficient=md.results.StressbalanceSolution.FrictionCoefficient;

	%Plot and save
	plotmodel(md,...
		'data',md.initialization.vel,'title','Observed velocity',...
		'data',md.results.StressbalanceSolution.Vel,'title','Modeled Velocity',...
		'data',md.geometry.base,'title','Bed elevation',...
		'data',md.results.StressbalanceSolution.FrictionCoefficient,'title','Friction Coefficient',...
		'colorbar#all','on','colorbartitle#1-2','(m/yr)',...
		'caxis#1-2',([1.5,4000]),...
		'colorbartitle#3','(m)', 'log#1-2',10);

	save ./Models/PIG_Control_drag md;
end 

if any(steps==6) %Transient Run #1 

	md = loadmodel('./Models/PIG_Control_drag');	

	md.inversion.iscontrol=0;
	md.transient.ismasstransport=1;
	md.transient.isstressbalance=1;
	md.transient.isgroundingline=1;
	md.transient.ismovingfront=0;
	md.transient.isthermal=0;
	md.verbose.solution=1;
	md.timestepping.time_step=0.1;
	md.timestepping.final_time=10;
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation'};

	%Set melt to 25 m/yr under floating ice
	md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
	md.basalforcings.floatingice_melting_rate=25*ones(md.mesh.numberofvertices,1);


	md=solve(md,'Transient');
	time      = cell2mat({md.results.TransientSolution(:).time});
	V_control = cell2mat({md.results.TransientSolution(:).IceVolumeAboveFloatation});
	plot(time,V_control); legend({'Control run'});

	% Save model
	save ./Models/PIG_Transient md;
end 

if any(steps==7) %High Melt #2
	md = loadmodel('./Models/PIG_Transient');	

	%Set melt to 60 m/yr under floating ice
	md.basalforcings.floatingice_melting_rate=60*ones(md.mesh.numberofvertices,1);

	md=solve(md,'Transient');
	V_melt = cell2mat({md.results.TransientSolution(:).IceVolumeAboveFloatation});
	plot(time,[V_control',V_melt']); legend({'Control run','High melt'});

	save ./Models/PIG_HighMelt md;
end 

if any(steps==8) %Ice Front retreat 
	md = loadmodel('./Models/PIG_Transient');	

	pos = find(ContourToNodes(md.mesh.x,md.mesh.y,'FrontRetreat.exp',2));
	md.mask.ice_levelset(pos) = +1; %Deactivate nodes in retreated area

	md=solve(md,'Transient');
	V_retreat = cell2mat({md.results.TransientSolution(:).IceVolumeAboveFloatation});
	plot(time,[V_control',V_melt',V_retreat']); legend({'Control run','High melt','Retreat'});

	save ./Models/PIG_FrontRetreat md;
end 

if any(steps==9) %High surface mass balance #3 
	%Load model from PIG_Transient
	%...

	%Change surface mass balance (x2)
	%...

	%Solve
	%...

	%Get volume time series
	%...

	%plot
	plot(time,[V_control',V_melt',V_retreat',V_smb']); legend({'Control run','High melt','Retreat','SMB'});

	%Save model
	save ./Models/PIG_SMB md;
end 
