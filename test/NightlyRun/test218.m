%Test Name: SquareShelfConstrainedDakotaB
md=squaremesh(model(),1000000,1000000,5,5);
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareShelf2.par');
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',3);

%redo the parameter file for this special shelf. 
%constant thickness, constrained (vy=0) flow into an icefront, 
%from 0 m/yr at the grounding line.

%needed later
ymin=min(md.mesh.y);
ymax=max(md.mesh.y);
xmin=min(md.mesh.x);
xmax=max(md.mesh.x);

di=md.materials.rho_ice/md.materials.rho_water;

h=1000;
md.geometry.thickness=h*ones(md.mesh.numberofvertices,1);
md.geometry.base=-md.materials.rho_ice/md.materials.rho_water*md.geometry.thickness;
md.geometry.surface=md.geometry.base+md.geometry.thickness;

%Initial velocity and pressure
md.initialization.vx=zeros(md.mesh.numberofvertices,1);
md.initialization.vy=zeros(md.mesh.numberofvertices,1);
md.initialization.vz=zeros(md.mesh.numberofvertices,1);
md.initialization.pressure=zeros(md.mesh.numberofvertices,1);

%Materials
md.initialization.temperature=(273-20)*ones(md.mesh.numberofvertices,1);
md.materials.rheology_B=paterson(md.initialization.temperature);
md.materials.rheology_n=3*ones(md.mesh.numberofelements,1);

%Boundary conditions:
md.stressbalance.spcvx=NaN*ones(md.mesh.numberofvertices,1);
md.stressbalance.spcvy=NaN*ones(md.mesh.numberofvertices,1);
md.stressbalance.spcvz=NaN*ones(md.mesh.numberofvertices,1);

%constrain flanks to 0 normal velocity
pos=find(md.mesh.x==xmin | md.mesh.x==xmax);
md.stressbalance.spcvx(pos)=0;
md.stressbalance.spcvz(pos)=NaN;

%constrain grounding line to 0 velocity
pos=find(md.mesh.y==ymin);
md.stressbalance.spcvx(pos)=0;
md.stressbalance.spcvy(pos)=0;

%partitioning
npart=md.mesh.numberofvertices;
partition=partitioner(md,'package','linear','npart',npart)-1;

%Dakota options

%dakota version
version=IssmConfig('_DAKOTA_VERSION_'); version=version(1:3); version=str2num(version);

%variables
md.qmu.variables.rheology_B=normal_uncertain('descriptor','scaled_MaterialsRheologyB',...
	'mean',ones(md.mesh.numberofvertices,1),...
	'stddev',.05*ones(md.mesh.numberofvertices,1),...
	'partition',partition);

%responses
md.qmu.responses.MaxVel=response_function('descriptor','MaxVel');

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
md.stressbalance.reltol=10^-10; %tighten for qmu analyses
md.qmu.isdakota=1;

%solve
md=solve(md,'Stressbalance','overwrite','y');

%Fields and tolerances to track changes
md.qmu.results=md.results.dakota;
md.results.dakota.importancefactors=importancefactors(md,'scaled_MaterialsRheologyB','MaxVel',partition)';
field_names     ={'importancefactors'};
field_tolerances={1e-10};
field_values={...
         md.results.dakota.importancefactors,...
	};
