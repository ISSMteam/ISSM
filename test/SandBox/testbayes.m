%Test Name: Test Queos
md=triangle(model(),'../Exp/Square.exp',150000.);
md=setmask(md,'../Exp/SquareShelf.exp','');
md=parameterize(md,'../Par/SquareSheetShelf.par');
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',3);

%constrain all velocities to 1 m/yr, in the y-direction
pos=find(md.mesh.y==0); 
md.stressbalance.spcvx(pos)=0;
md.stressbalance.spcvy(pos)=0;
%md.stressbalance.spcvy(:)=1;
md.stressbalance.spcvz(:)=0;

%Dakota options
md.qmu.isdakota=1;

%partitioning
md.qmu.numberofpartitions=md.mesh.numberofvertices;
md=partitioner(md,'package','linear','npart',md.qmu.numberofpartitions);
md.qmu.vpartition=md.qmu.vpartition-1;


%input/output:
md.qmu.variables.drag_coefficient=normal_uncertain('scaled_FrictionCoefficient',1,0.10);
md.qmu.responses.misfit=calibration_function('MaxVel');
%md.qmu.responses.MaxVel=response_function('MaxVel',[],[0.0001 0.001 0.01 0.25 0.5 0.75 0.99 0.999 0.9999]);

%bayesian calibration:
md.qmu.method     =dakota_method('bayes_calibration');
md.qmu.method(end)=dmeth_params_set(md.qmu.method(end),...
'seed',1234,...
'samples',2000,...
'queso',true,...
'metropolis_hastings',true,...
'proposal_covariance',true,...
'diagonal',true,...
'values',[num2str(numel(fields(md.qmu.variables))*md.qmu.numberofpartitions) ' * 0.00025']);

%direct solve: 
md.qmu.params.direct=true;
md.qmu.params.interval_type='forward';
md.qmu.params.tabular_graphics_data=true;
md.qmu.params.analysis_driver='matlab';
md.qmu.params.evaluation_scheduling='master';
md.qmu.params.processors_per_evaluation=2;

md.qmu.isdakota=1;

md=solve(md,'Stressbalance','overwrite','y');

%Fields and tolerances to track changes
md.qmu.results=md.results.dakota;
