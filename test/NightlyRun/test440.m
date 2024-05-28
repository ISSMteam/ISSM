%Test Name: SquareSheetShelfDakotaScaledResponseLinearPart
md=triangle(model(),'../Exp/Square.exp',200000.);
md=setmask(md,'../Exp/SquareShelf.exp','');
md=parameterize(md,'../Par/SquareSheetShelf.par');
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',3);

%partitioning
npart = md.mesh.numberofvertices;
partition=partitioner(md,'package','linear','npart',npart)-1;
md.qmu.isdakota=1;

%Dakota options

%dakota version
version=IssmConfig('_DAKOTA_VERSION_'); version=version(1:3); version=str2num(version);

%variables
md.qmu.variables.rho_ice=normal_uncertain('descriptor','MaterialsRhoIce','mean',1,'stddev',.01);

%responses
md.qmu.responses.MaxVel=response_function('descriptor','scaled_Thickness','partition',partition);

%method
md.qmu.method     =dakota_method('nond_l');

%parameters
md.qmu.params.direct=true;
md.qmu.params.interval_type='forward';

if version>=6,
	md.qmu.params.analysis_driver='matlab';
	md.qmu.params.evaluation_scheduling='master';
	md.qmu.params.processors_per_evaluation=2;
else
	md.qmu.params.analysis_driver='stressbalance';
	md.qmu.params.evaluation_concurrency=1;
end

%imperative! 
md.stressbalance.reltol=10^-5; %tighten for qmu analysese

%solve
md=solve(md,'Stressbalance','overwrite','y');
md.qmu.results=md.results.dakota;

%test on thickness
h=zeros(npart,1);
for i=1:npart,
	h(i)=md.qmu.results.dresp_out(i).mean;
end

%project onto grid
thickness=h(partition+1);

%Fields and tolerances to track changes
field_names     ={'Thickness'};
field_tolerances={1e-10};
field_values={thickness};
