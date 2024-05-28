steps=[1];

if any(steps==1) %Mesh Generation #1 
	%Mesh parameters
	domain =['./DomainOutline.exp'];
	hinit=10000;	% element size for the initial mesh
	hmax=40000;		% maximum element size of the final mesh
	hmin=5000;		% minimum element size of the final mesh
	gradation=1.7;	% maximum size ratio between two neighboring elements
	err=8;			% maximum error between interpolated and control field

	% Generate an initial uniform mesh (resolution = hinit m)
	md=bamg(model,'domain',domain,'hmax',hinit);

	% Load Velocities
	nsidc_vel='../Data/Antarctica_ice_velocity.nc';

	% Get necessary data to build up the velocity grid
	xmin	= ncreadatt(nsidc_vel,'/','xmin');
	ymax	= ncreadatt(nsidc_vel,'/','ymax');
	spacing	= ncreadatt(nsidc_vel,'/','spacing');
	nx		= double(ncreadatt(nsidc_vel,'/','nx'));
	ny		= double(ncreadatt(nsidc_vel,'/','ny'));
	vx		= double(ncread(nsidc_vel,'vx'));
	vy		= double(ncread(nsidc_vel,'vy'));

	% Read coordinates
	xmin = strtrim(xmin);
	xmin = str2num(xmin(1:end-2));
	ymax = strtrim(ymax);
	ymax = str2num(ymax(1:end-2));
	spacing = strtrim(spacing);
	spacing = str2num(spacing(1:end-2));

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

	%ploting
	plotmodel(md,'data','mesh')

	% Save model
	save ./Models/PIG_Mesh_generation md;
end 

if any(steps==2) %Masks #2 
	md = loadmodel('./Models/PIG_Mesh_generation');

	% Load SeaRISe dataset for Antarctica
	% http://websrv.cs.umt.edu/isis/index.php/Present_Day_Antarctica
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

	%ploting
	plotmodel(md,'data',md.mask.ocean_levelset,'title','grounded/floating','data',md.mask.ice_levelset,'title','ice/no-ice')

	% Save model
	save ./Models/PIG_SetMask md;
end 

if any(steps==3) %Parameterization #3 
	md = loadmodel('./Models/PIG_SetMask');
	md = parameterize(md,'./Pig.par');

	% Use a MacAyeal flow model
	md = setflowequation(md,'SSA','all');

	% Save model
	save ./Models/PIG_Parameterization md;
end 

if any(steps==4) %Control Method #4 
	md = loadmodel('./Models/PIG_Parameterization');

	% Control general
	md.inversion.iscontrol=1;
	md.inversion.maxsteps=20;
	md.inversion.maxiter=40;
	md.inversion.dxmin=0.1;
	md.inversion.gttol=1.0e-4;
	md.verbose=verbose('control',true);

	% Cost functions
	md.inversion.cost_functions=[101 103 501];
	md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,3);
	md.inversion.cost_functions_coefficients(:,1)=1;
	md.inversion.cost_functions_coefficients(:,2)=1;
	md.inversion.cost_functions_coefficients(:,3)=8e-15;

	% Controls
	md.inversion.control_parameters={'FrictionCoefficient'};
	md.inversion.min_parameters=1*ones(md.mesh.numberofvertices,1);
	md.inversion.max_parameters=200*ones(md.mesh.numberofvertices,1);

	% Additional parameters
	md.stressbalance.restol=0.01;
	md.stressbalance.reltol=0.1;
	md.stressbalance.abstol=NaN;

	% Solve
	md.toolkits=toolkits;
	md.cluster=generic('name',oshostname,'np',2);
	md=solve(md,'Stressbalance');

	% Update model friction fields accordingly
	md.friction.coefficient=md.results.StressbalanceSolution.FrictionCoefficient;

	plotmodel(md,'data',md.friction.coefficient)

	% Save model
	save ./Models/PIG_Control_drag md;
end 

if any(steps==5) %Plot #5 
	md = loadmodel('./Models/PIG_Control_drag');

	plotmodel(md,...
		'data',md.initialization.vel,'title','Observed velocity',...
		'data',md.results.StressbalanceSolution.Vel,'title','Modeled Velocity',...
		'data',md.geometry.base,'title','Bed elevation',...
		'data',md.results.StressbalanceSolution.FrictionCoefficient,'title','Friction Coefficient',...
		'colorbar#all','on','colorbartitle#1-2','(m/yr)',...
		'caxis#1-2',([1.5,4000]),...
		'colorbartitle#3','(m)', 'log#1-2',10);
end 

if any(steps==6) %Higher-Order #6 
	% Load Model

	% Disable inversion

	% Extrude Mesh

	% Set Flowequation

	% Solve

	% Save Model

end % step 6 end

if any(steps==7) %Plot #7 
	mdHO = loadmodel('./Models/PIG_ModelHO');
	mdSSA = loadmodel('./Models/PIG_Control_drag');

	basal=find(mdHO.mesh.vertexonbase);
	surf=find(mdHO.mesh.vertexonsurface);

	plotmodel(mdHO,'nlines',3,'ncols',2,'axis#all','equal',...
		'data',mdHO.initialization.vel,'title','Observed velocity',...
		'data',(mdHO.results.StressbalanceSolution.Vel(surf)-mdHO.initialization.vel(surf)),'title','(HO-observed) velocities',...
		'data',mdSSA.results.StressbalanceSolution.Vel,'title','Modeled SSA Velocity',...
		'data',(mdHO.results.StressbalanceSolution.Vel(surf)-mdSSA.results.StressbalanceSolution.Vel),'title','(HO-SSA) velocities',...
		'data',mdHO.results.StressbalanceSolution.Vel,'title','Modeled HO surface Velocities',...
		'data',(mdHO.results.StressbalanceSolution.Vel(surf)-mdHO.results.StressbalanceSolution.Vel(basal)),'title','(HOsurf-HO base) velocities',...
		'caxis#1',([1.5,4000]),'caxis#3',([1.5,4000]),'caxis#5',([1.5,4000]),...
		'colorbar#all','on','view#all',2,...
		'colorbartitle#all','(m/yr)',...
		'layer#5',1, 'log#1', 10,'log#3', 10,'log#5', 10);
end 
