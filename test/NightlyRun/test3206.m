%Test Name: SquareShelfTransientCalibrationWithParamcodipack

%Generate observations
md = model;
md=triangle(model(),'../Exp/Square.exp',50000.);
md = setmask(md,'all','');
md = parameterize(md,'../Par/SquareShelf.par');
md = setflowequation(md,'SSA','all');
md.cluster = generic('np',2);

%Create real time series for B
md.timestepping.interp_forcing = 0;
md.timestepping.final_time = 2*md.timestepping.time_step;
md.materials.rheology_B = 1.8e8*ones(md.mesh.numberofvertices,2);
md.materials.rheology_B(find(md.mesh.x<md.mesh.y),2)=1.4e8;
md.materials.rheology_B=[md.materials.rheology_B;0.01 2*md.timestepping.time_step];

%Initial values
md.initialization.vx = zeros(md.mesh.numberofvertices,1);
md.initialization.vy = zeros(md.mesh.numberofvertices,1);
md.initialization.pressure = zeros(md.mesh.numberofvertices,1);
md.initialization.temperature = zeros(md.mesh.numberofvertices,1);
md.basalforcings.geothermalflux = zeros(md.mesh.numberofvertices,1);
md.thermal.spctemperature = NaN(md.mesh.numberofvertices,1);

%Param
md.basalforcings=linearbasalforcings();
md.basalforcings.deepwater_melting_rate=50.; % m/yr ice equivalent
md.basalforcings.deepwater_elevation=-500;
md.basalforcings.upperwater_melting_rate=0; % no melting for zb>=0
md.basalforcings.upperwater_elevation=0; % sea level
md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1); % no melting on grounded ice
md.basalforcings.perturbation_melting_rate(:)=0;
md.transient.isthermal = 0;

md = solve(md,'tr');

%Set cost function
count = 1;
for i=1:numel(md.results.TransientSolution)
	vx_obs = md.results.TransientSolution(i).Vx;
	vy_obs = md.results.TransientSolution(i).Vy;
	time   = md.results.TransientSolution(i).time;
	weights= ones(md.mesh.numberofvertices,1);

	md.outputdefinition.definitions{count}=cfsurfacelogvel('name',['LogVelMis' num2str(count)],...
		'definitionstring',['Outputdefinition' num2str(count)],...
		'vxobs_string','VxObs','vxobs',vx_obs,...
		'vyobs_string','VyObs','vyobs',vy_obs,...
		'weights',weights,'weights_string','WeightsSurfaceObservation',...
		'datatime',time);
	md.autodiff.dependents{count} = dependent('name',['Outputdefinition' num2str(count)],'type','scalar','fos_reverse_index',1);
	count = count+1;
end

%Deal with vx separately
vx_obs  = [[md.results.TransientSolution(:).Vx]/md.constants.yts; [md.results.TransientSolution(:).time]];
weights = [ones(size(vx_obs,1)-1,1); 0];
md.outputdefinition.definitions{count}=cfsurfacesquaretransient('name',['VxMisfit_Transient'],...
	'definitionstring',['Outputdefinition' num2str(count)],...
	'model_string','Vx','observations_string','VxObs',...
	'observations',vx_obs,'weights',500*weights,'weights_string','WeightsSurfaceObservation');
md.autodiff.dependents{count} = dependent('name',['Outputdefinition' num2str(count)],'type','scalar','fos_reverse_index',1);
count = count+1;

vy_obs  = [[md.results.TransientSolution(:).Vy]/md.constants.yts; [md.results.TransientSolution(:).time]];
md.outputdefinition.definitions{count}=cfsurfacesquaretransient('name',['VyMisfit_Transient'],...
	'definitionstring',['Outputdefinition' num2str(count)],...
	'model_string','Vy','observations_string','VyObs',...
	'observations',vy_obs,'weights',weights,'weights_string','WeightsSurfaceObservation');
md.autodiff.dependents{count} = dependent('name',['Outputdefinition' num2str(count)],'type','scalar','fos_reverse_index',1);
count = count+1;

surf_obs  = [[md.results.TransientSolution(:).Surface]; [md.results.TransientSolution(:).time]];
md.outputdefinition.definitions{count}=cfsurfacesquaretransient('name',['SurfMisfit_Transient'],...
	'definitionstring',['Outputdefinition' num2str(count)],...
	'model_string','Surface','observations_string','SurfaceObservation',...
	'observations',surf_obs,'weights',weights/(md.constants.yts),'weights_string','WeightsSurfaceObservation');
md.autodiff.dependents{count} = dependent('name',['Outputdefinition' num2str(count)],'type','scalar','fos_reverse_index',1);
count = count+1;

%Independent
md.materials.rheology_B(1:end-1,:) = 1.8e8;
min_params = md.materials.rheology_B; min_params(1:end-1,:) = cuffey(273);
max_params = md.materials.rheology_B; max_params(1:end-1,:) = cuffey(200);
md.autodiff.independents{1} = independent('name','MaterialsRheologyBbar',...
	'control_size',size(md.materials.rheology_B,2),...
	'type','vertex',... %Really needed??
	'min_parameters',min_params,...
	'max_parameters',max_params,...
	'control_scaling_factor',1e8);

md.basalforcings.deepwater_melting_rate=1.; % m/yr ice equivalent
field =md.basalforcings.deepwater_melting_rate/md.constants.yts;
name = 'BasalforcingsDeepwaterMeltingRate';
scaling = 50/md.constants.yts;
md.autodiff.independents{2} = independent('name',name,'type','vertex','nods',md.mesh.numberofvertices,...
	'control_size', size(field,2), 'min_parameters',1e-5*field, 'max_parameters',100*field, 'control_scaling_factor',scaling);

md.inversion=adm1qn3inversion(md.inversion);
md.inversion.iscontrol=1;
md.inversion.maxiter=3;
md.inversion.maxsteps=md.inversion.maxiter;
md.inversion.dxmin=1e-5;
md.autodiff.isautodiff=1;
md.autodiff.driver='fos_reverse';
md.settings.checkpoint_frequency = 2;

%Go solve!
md.verbose=verbose(0);
md=solve(md,'tr');

%Fields and tolerances to track changes
field_names     ={'Gradient1','Gradient2','Misfit','Rheology','DeepMelt'};
field_tolerances={1e-10,1e-10,1e-10,1e-10,1e-10};
field_values={...
	(md.results.TransientSolution(1).Gradient1),...
	(md.results.TransientSolution(1).Gradient2),...
	(md.results.TransientSolution(1).J),...
	(md.results.TransientSolution(1).MaterialsRheologyBbar),...
	(md.results.TransientSolution(1).BasalforcingsDeepwaterMeltingRate),...
	};


return;
%The code below validates the gradient, run only with maxiter=1 above!
disp('Testing Gradient');
index = 3;
dJdB_ad = md.results.TransientSolution(1).Gradient1(index);
delta=0.001;
B1=md.materials.rheology_B(index);
%B1=md.basalforcings.deepwater_melting_rate;
B0=B1*(1.-delta);
B2=B1*(1.+delta);
deltaB=(B2-B0);

list = {}; for i=1:numel(md.outputdefinition.definitions), list{i} = md.autodiff.dependents{i}.name;end
md.transient.requested_outputs = list;
md.autodiff.isautodiff=false;
md.inversion.iscontrol=false;
md2=md;

%forward
md=md2;
md.materials.rheology_B(index)=B0;
%md.basalforcings.deepwater_melting_rate = B0;
md=solve(md,'tr');
J0 = 0;
for i=1:numel(md.outputdefinition.definitions), eval(['J0 = J0 + md.results.TransientSolution(end).' list{i} ';']); end

%backward
md=md2;
md.materials.rheology_B(index)=B2;
%md.basalforcings.deepwater_melting_rate = B2;
md=solve(md,'tr');
J2 = 0;
for i=1:numel(md.outputdefinition.definitions), eval(['J2 = J2 + md.results.TransientSolution(end).' list{i} ';']); end

%compute resulting derivative
dJdB_an=(J2-J0)/deltaB;

disp(sprintf('dJ/dB: analytical:  %16.16g\n       using ad:    %16.16g\n',dJdB_an,dJdB_ad));
