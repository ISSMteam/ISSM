%Test Name: SquareShelfTransientCalibrationNBVcodipack

%Generate observations
md = model;
md=triangle(model(),'../Exp/Square.exp',100000.);
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

md = solve(md,'tr');
%plotmodel(md,'axis#all','tight','data',md.materials.rheology_B(1:end-1,1),'caxis#all',[ 1.3 1.9]*10^8,'title','"True" B',...
%'data',md.materials.rheology_B(1:end-1,2),'title','"True" B 2')

%Modify rheology, now constant
md.materials.rheology_B(1:end-1,:) = 1.8e8;

%Set cost function
count = 1;
for i=1:numel(md.results.TransientSolution)
	vx_obs = md.results.TransientSolution(i).Vx;
	vy_obs = md.results.TransientSolution(i).Vy;
	z_obs  = md.results.TransientSolution(i).Surface;

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

	md.outputdefinition.definitions{count}=cfsurfacesquare('name',['VyMisfit' num2str(count)],...
		'definitionstring',['Outputdefinition' num2str(count)],...
		'model_string','Vy','observation_string','VyObs',...
		'observation',vy_obs/md.constants.yts,'weights',weights,'weights_string','WeightsSurfaceObservation',...
		'datatime',time);
	md.autodiff.dependents{count} = dependent('name',['Outputdefinition' num2str(count)],'type','scalar','fos_reverse_index',1);
	count = count+1;

	md.outputdefinition.definitions{count}=cfsurfacesquare('name',['VxMisfit' num2str(count)],...
		'definitionstring',['Outputdefinition' num2str(count)],...
		'model_string','Vx','observation_string','VxObs',...
		'observation',vx_obs/md.constants.yts,'weights',500*weights,'weights_string','WeightsSurfaceObservation',...
		'datatime',time);
	md.autodiff.dependents{count} = dependent('name',['Outputdefinition' num2str(count)],'type','scalar','fos_reverse_index',1);
	count = count+1;

	md.outputdefinition.definitions{count}=cfsurfacesquare('name',['DEMMisfit' num2str(count)],...
		'definitionstring',['Outputdefinition' num2str(count)],...
		'model_string','Surface','observation_string','SurfaceObservation',...
		'observation',z_obs,...
		'weights',1/(md.constants.yts)*weights,...
		'weights_string','WeightsSurfaceObservation',...
		'datatime',time);
	md.autodiff.dependents{count} = dependent('name',['Outputdefinition' num2str(count)],'type','scalar','fos_reverse_index',1);
	count = count+1;
end

%Independent
min_params = md.materials.rheology_B; min_params(1:end-1,:) = cuffey(273);
max_params = md.materials.rheology_B; max_params(1:end-1,:) = cuffey(200);
md.autodiff.independents{1} = independent('name','MaterialsRheologyBbar',...
	'control_size',size(md.materials.rheology_B,2),...
	'type','vertex',... %Really needed??
	'min_parameters',min_params,...
	'max_parameters',max_params,...
	'control_scaling_factor',1e8);

md.inversion=adm1qn3inversion(md.inversion);
md.inversion.iscontrol=1;
md.inversion.maxiter=4;
md.inversion.maxsteps=md.inversion.maxiter;
md.inversion.dxmin=1e-5;
md.autodiff.isautodiff=1;
md.autodiff.driver='fos_reverse';

%Go solve!
md.verbose=verbose(0);
md=solve(md,'tr');
%plotmodel(md,'axis#all','tight','data',md.results.TransientSolution(1).MaterialsRheologyBbar(:,1),'caxis#all',[ 1.3 1.9]*10^8,'title','B1',...
%'data',md.results.TransientSolution(1).MaterialsRheologyBbar(:,2),'title','B2')

%Fields and tolerances to track changes
field_names     ={'Gradient','Misfit','Rheology'};
field_tolerances={2e-12,1e-12,1e-12};
field_values={...
	(md.results.TransientSolution(1).Gradient1),...
	(md.results.TransientSolution(1).J),...
	(md.results.TransientSolution(1).MaterialsRheologyBbar),...
	};
