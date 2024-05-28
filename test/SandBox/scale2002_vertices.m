%Test Name: SquareSheetShelfDiadSSA3dDakota
md=triangle(model(),'../Exp/Square.exp',300000.);
md=setmask(md,'../Exp/SquareShelf.exp','');
md=parameterize(md,'../Par/SquareSheetShelf.par');
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',2);

%partitioning
partition=partitioner(md,'package','linear','npart',md.mesh.numberofvertices)-1;
md.qmu.isdakota=1;

nsteps=100;
md.friction.coefficient=ones(md.mesh.numberofvertices+1,nsteps);
md.friction.coefficient(end,:)=2000.5:1:2099.5;

%variables
md.qmu.variables.drag_coefficient=normal_uncertain('descriptor','scaled_FrictionCoefficient',...
	'mean',ones(md.mesh.numberofvertices,nsteps),...
	'stddev',.01*ones(md.mesh.numberofvertices,nsteps),...
	'partition',partition,'nsteps',nsteps);

%responses
md.qmu.responses.MaxVel=response_function('descriptor','MaxVel');

%algorithm: 
md.qmu.method     =dakota_method('nond_samp');
md.qmu.method(end)=dmeth_params_set(md.qmu.method(end),...
	'seed',1234,...
	'samples',3,...
	'sample_type','random');

%parameters
md.qmu.params.direct=true;
md.qmu.params.interval_type='forward';
md.qmu.params.analysis_driver='matlab';
md.qmu.params.evaluation_scheduling='master';
md.qmu.params.processors_per_evaluation=1;
md.qmu.params.tabular_graphics_data=true;
md.qmu.output=1;

md.verbose=verbose(0);

md.timestepping.start_time=1999.5;
md.timestepping.interp_forcings=0;
md.timestepping.time_step=1;
md.timestepping.final_time=2099.5
md.transient.issmb=0;
md.transient.ismasstransport=0;
md.transient.isthermal=0;

md.stressbalance.requested_outputs{end+1}='FrictionCoefficient';

%imperative! 
md.stressbalance.reltol=10^-5; %tighten for qmu analysese

%solve
md=solve(md,'tr','overwrite','y');

%check scaling worked ok: 
n=md.qmu.method.params.samples;
field_values={};

for i=1:n,
	md2=model(); md2.results=md.results.dakota.modelresults{i};
	uq=md2.results.TransientSolution(1).uq_variables;
	
	h2=ones(md.mesh.numberofvertices+1,100);

	for j=1:md.mesh.numberofvertices,
		partj=uq((j-1)*100+1 : j*100)
		h2(j,:)=h2(j,:)'.*partj;
	end
	
	h=resultstomatrix(md2,'TransientSolution','FrictionCoefficient');
	h2(end,:)=h(end,:);
	size(find(h2-h))

	field_values{end+1}=size(find(h2-h));
end


