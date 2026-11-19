function result = runme(stage,np,seed)% {{{
%RUNME Train the native friction network through an ISSM solve.
%   runme('stressbalance') checks the gradient and fits one SSA solve.
%   runme('transient') checks the gradient and fits a two-step rollout.
%   An optional second argument sets the MPI rank count (default 1).
%   An optional third argument sets the initialization seed (default 1).

	%% Options and paths
	if nargin<1, stage='stressbalance'; end
	if nargin<2, np=1; end
	if nargin<3, seed=1; end
	stage=validatestring(stage,{'stressbalance','transient'});
	example_dir=fileparts(mfilename('fullpath'));
	original_dir=pwd;
	restore_dir=onCleanup(@() cd(original_dir));
	cd(fullfile(example_dir,'..','..','test','NightlyRun'));

	%% Build the grounded SSA model
	C=100;
	md=triangle(model(),'../Exp/Square.exp',200000);
	md=setmask(md,'','');
	md=parameterize(md,'../Par/SquareSheetConstrained.par');
	md=setflowequation(md,'SSA','all');
	md.cluster=generic('name',oshostname(),'np',np);
	md.stressbalance.restol=1e-4;
	md.stressbalance.reltol=1e-2;
	md.stressbalance.abstol=10;
	md.stressbalance.maxiter=50;
	md.settings.checkpoint_frequency=0;
	md.friction=frictionweertman();
	md.friction.C=C*ones(md.mesh.numberofvertices,1);
	md.friction.m=3*ones(md.mesh.numberofelements,1);

	%% Configure the solve and generate Weertman reference data
	switch stage
		case 'stressbalance'
			solution='Stressbalance';
		case 'transient'
			solution='Transient';
			md.transient=deactivateall(md.transient);
			md.transient.issmb=1;
			md.transient.isstressbalance=1;
			md.transient.ismasstransport=1;
			md.timestepping.time_step=0.05;
			md.timestepping.final_time=md.timestepping.start_time+2*md.timestepping.time_step;
			md.settings.output_frequency=1;
	end
	result_field=[solution 'Solution'];
	prefix=['NN' solution];
	truth=run_solve(md,solution,[prefix 'Truth']);
	reference=truth.results.(result_field);
	assert(numel(reference)==1+strcmp(stage,'transient'),'Incomplete reference rollout.');
	assert_converged(reference,md.stressbalance.maxiter);
	speed=[reference.Vel];
	assert(all(isfinite(speed(:))) && max(speed(:))>0,'Invalid reference speeds.');

	%% Initialize and run the neural friction model
	md.friction=frictionnn();
	md.friction.C=truth.friction.C;
	theta=initial_network(seed);
	normalization=network_normalization(C,speed,[reference.FrictionAlpha2],md.constants.yts);
	md.friction.nn_weights=theta;
	md.friction.nn_normalization=normalization;
	initial=run_solve(md,solution,[prefix 'Initial']);
	initial_solution=initial.results.(result_field);
	assert_converged(initial_solution,md.stressbalance.maxiter);

	%% Define the objective and configure CoDiPack
	% Use the same initial state and fixed observations for every evaluation.
	area=sum(GetAreas(md.mesh.elements,md.mesh.x,md.mesh.y));
	speed_scale=max(max(speed(:)),1)/md.constants.yts;
	weights=ones(md.mesh.numberofvertices,1)/(area*speed_scale^2*numel(reference));
	for k=1:numel(reference)
		definition=sprintf('Outputdefinition%d',k);
		time=md.timestepping.start_time;
		if strcmp(stage,'transient'), time=reference(k).time; end
		md.outputdefinition.definitions{k}=cfsurfacesquare( ...
			'name',sprintf('SpeedMisfit%d',k),'definitionstring',definition, ...
			'model_string','Vel','observation_string','SurfaceObservation', ...
			'observation',reference(k).Vel/md.constants.yts,'weights',weights,'datatime',time);
		md.autodiff.dependents{k}=dependent('name',definition,'type','scalar','fos_reverse_index',1);
	end
	md.autodiff.independents={independent('name','FrictionNNWeights','type','scalar', ...
		'control_size',numel(theta),'min_parameters',theta-10,'max_parameters',theta+10, ...
		'control_scaling_factor',1)};
	md.autodiff.isautodiff=1;
	md.autodiff.driver='fos_reverse';
	md.inversion=adm1qn3inversion();
	md.inversion.iscontrol=1;
	md.inversion.maxsteps=100;
	md.inversion.maxiter=300;
	md.inversion.dxmin=1e-8;
	md.inversion.dfmin_frac=0.1;
	md.inversion.gttol=1e-6;

	%% Check the gradient against finite differences
	% A tolerance of 1 lets M1QN3 return after the initial AD evaluation.
	probe=md;
	probe.inversion.gttol=1;
	base=evaluate(probe,solution,result_field,[prefix 'Gradient']);
	gradient=base.Gradient1;
	% Check the final output bias and four other parameters with large gradients.
	[~,ranked]=sort(abs(gradient(1:end-1)),'descend');
	indices=[numel(theta),reshape(ranked(1:4),1,[])];
	errors=zeros(size(indices));
	fd=zeros(size(indices));
	h=1e-5;
	for j=1:numel(indices)
		k=indices(j);
		delta=zeros(size(theta));
		delta(k)=h;
		probe.friction.nn_weights=theta+delta;
		plus=evaluate(probe,solution,result_field,sprintf('%sPlus%d',prefix,k));
		probe.friction.nn_weights=theta-delta;
		minus=evaluate(probe,solution,result_field,sprintf('%sMinus%d',prefix,k));
		fd(j)=(plus.J(1,end)-minus.J(1,end))/(2*h);
		ad=gradient(k);
		errors(j)=abs(fd(j)-ad)/max([abs(fd(j)),abs(ad),1e-12]);
	end

	for j=1:numel(indices)
		k=indices(j);
		fprintf('%s gradient check, parameter %d: AD %.8g, FD %.8g, relative error %.3g\n', ...
			stage,k,gradient(k),fd(j),errors(j));
	end
	gradient_error=max(errors);
	assert(all(isfinite(errors)) && abs(gradient(end))>1e-12 && gradient_error<1e-2, ...
		'Neural friction gradient check failed.');

	%% Fit the neural friction weights
	fit=run_solve(md,solution,[prefix 'Fit']);
	fitted=fit.results.(result_field);
	trained_weights=reshape(fitted(1).FrictionNNWeights,1,[]);
	objective=fitted(1).J(:,end);
	assert(all(isfinite(trained_weights)) && min(objective)<objective(1), ...
		'Neural friction training did not reduce the objective.');

	%% Validate the fitted weights in a forward run
	% Exercise the forward law using only the exported weights.
	md.friction.nn_weights=trained_weights;
	md.inversion.iscontrol=0;
	md.autodiff.isautodiff=0;
	trained=run_solve(md,solution,[prefix 'Trained']);
	prediction=trained.results.(result_field);
	assert_converged(prediction,md.stressbalance.maxiter);
	assert(numel(prediction)==numel(reference),'Incomplete forward rollout.');
	initial_rmse=sqrt(mean(([initial_solution.Vel]-speed).^2,'all'));
	trained_rmse=sqrt(mean(([prediction.Vel]-speed).^2,'all'));
	fprintf('%s speed RMSE: %.6g -> %.6g m/yr\n',stage,initial_rmse,trained_rmse);
	fprintf('%s objective: %.8g -> %.8g\n',stage,objective(1),min(objective));
	initial_drag_error=drag_error(initial_solution,C,md.constants.yts);
	trained_drag_error=drag_error(prediction,C,md.constants.yts);
	fprintf('%s drag relative L2 error at sampled speeds: %.4g -> %.4g\n', ...
		stage,initial_drag_error,trained_drag_error);

	%% Save the results
	result=struct('seed',seed,'trained_weights',trained_weights,'normalization',normalization, ...
		'objective',objective,'gradient_error',gradient_error, ...
		'initial_rmse',initial_rmse,'trained_rmse',trained_rmse, ...
		'initial_drag_error',initial_drag_error,'trained_drag_error',trained_drag_error, ...
		'reference',reference,'prediction',prediction);
	output_dir=fullfile(example_dir,'Results');
	if ~exist(output_dir,'dir'), mkdir(output_dir); end
	save(fullfile(output_dir,sprintf('nn_%s_np%d.mat',stage,np)),'result');
end% }}}

function md=run_solve(md,solution,name)% {{{
	md.miscellaneous.name=name;
	if ~md.inversion.iscontrol
		if strcmp(solution,'Transient')
			md.transient.requested_outputs{end+1}='FrictionAlpha2';
		else
			md.stressbalance.requested_outputs{end+1}='FrictionAlpha2';
		end
	end
	md=solve(md,solution);
end% }}}
function first=evaluate(md,solution,result_field,name)% {{{
	md=run_solve(md,solution,name);
	results=md.results.(result_field);
	first=results(1);
end% }}}
function assert_converged(results,limit)% {{{
	assert(isfield(results,'StressbalanceConvergenceNumSteps'), ...
		'Stress balance iteration count was not returned.');
	iterations=[results.StressbalanceConvergenceNumSteps];
	% Transient iteration counts can be written without a step number. The
	% solver also reports any step that hits the limit in the output log.
	assert(~isempty(iterations) && all(iterations<limit) && ...
		~any(contains(string(results(1).outlog),'maximum number of nonlinear iterations')), ...
		'Stress balance reached its iteration limit.');
end% }}}
function err=drag_error(results,C,yts)% {{{
	% Evaluate the analytical law at the emulator run's own basal speeds.
	speed=[results.Vel]/yts;
	alpha2=[results.FrictionAlpha2];
	positive=speed>0;
	exact=C^2*speed(positive).^(-2/3);
	err=norm(alpha2(positive)-exact)/norm(exact);
end% }}}
function normalization=network_normalization(C,speed,alpha2,yts)% {{{
	% Use reference speeds and drag to set the input and output scales.
	positive_speed=speed(speed>0)/yts;
	speed_floor=max(min(positive_speed)/2,0.01/yts)/10;
	root_speed=sqrt(sqrt(positive_speed.^2+speed_floor^2));
	root_drag=sqrt(alpha2(isfinite(alpha2) & alpha2>0));
	assert(~isempty(root_drag),'No positive reference drag for normalization.');
	drag_scale=median(root_drag);
	normalization=[C mean(root_speed) C max(std(root_speed),eps) ...
		drag_scale drag_scale speed_floor];

end% }}}
function theta=initial_network(seed)% {{{
	stream=RandStream('mt19937ar','Seed',seed);
	layer_sizes=[2 32 32 1];
	theta=[];
	for k=1:numel(layer_sizes)-1
		nin=layer_sizes(k);
		nout=layer_sizes(k+1);
		limit=sqrt(6/(nin+nout));
		weights=limit*(2*rand(stream,nout,nin)-1);
		layer=[weights,zeros(nout,1)];
		theta=[theta,reshape(layer.',1,[])];
	end
end% }}}
