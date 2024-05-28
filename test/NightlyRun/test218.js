//Test Name: SquareShelfConstrainedDakotaB
var md = new model();
function zeros(...args) {
	var array = [];
	for (var i = 0; i < args[0]; ++i) {
		array.push(args.length == 1 ? 0 : zeros(args.slice(1)));
	}
	return array;
}
var md = new model();
squaremesh(model(),1000000,1000000,5,5);
setmask(md,'all','');
parameterize(md);
setflowequation(md,'SSA','all');
//md.cluster=generic('name',oshostname(),'np',3);

//redo the parameter file for this special shelf. 
//constant thickness, constrained (vy=0) flow into an icefront, 
//from 0 m/yr at the grounding line.

//needed later
ymin=Math.min.apply(null, md.mesh.y);
ymax=Math.max.apply(null, md.mesh.y);
xmin=Math.min.apply(null, md.mesh.x);
xmax=Math.max.apply(null, md.mesh.x);

di=md.materials.rho_ice/md.materials.rho_water;

h=1000;
md.geometry.thickness=h*ones(md.mesh.numberofvertices,1);
md.geometry.base=-md.materials.rho_ice/md.materials.rho_water*md.geometry.thickness;
md.geometry.surface=md.geometry.base+md.geometry.thickness;

//Initial velocity and pressure
md.initialization.vx=zeros(md.mesh.numberofvertices,1);
md.initialization.vy=zeros(md.mesh.numberofvertices,1);
md.initialization.vz=zeros(md.mesh.numberofvertices,1);
md.initialization.pressure=zeros(md.mesh.numberofvertices,1);

//Materials
md.initialization.temperature=(273-20)*ones(md.mesh.numberofvertices,1);
md.materials.rheology_B=paterson(md.initialization.temperature);
md.materials.rheology_n=3*ones(md.mesh.numberofelements,1);

//Boundary conditions:
md.stressbalance.spcvx=null*ones(md.mesh.numberofvertices,1);
md.stressbalance.spcvy=null*ones(md.mesh.numberofvertices,1);
md.stressbalance.spcvz=null*ones(md.mesh.numberofvertices,1);

//constrain flanks to 0 normal velocity
pos=[];
for (var i = 0; i < md.mesh.x==xmin | md.mesh.x==xmax.length; ++i) {
if ((md.mesh.x==xmin | md.mesh.x==xmax[i] !== 0)) {
		pos.push(i);
};
md.stressbalance.spcvx(pos)=0;
md.stressbalance.spcvz(pos)=null;

//constrain grounding line to 0 velocity
pos=[];
for (var i = 0; i < md.mesh.y==ymin.length; ++i) {
if ((md.mesh.y==ymin[i] !== 0)) {
		pos.push(i);
};
md.stressbalance.spcvx(pos)=0;
md.stressbalance.spcvy(pos)=0;

//partitioning
npart=md.mesh.numberofvertices;
partitioner(md,'package','linear','npart',npart);
md.qmu.partition=md.qmu.partition-1;

//Dakota options

//dakota version
version=IssmConfig('_DAKOTA_VERSION_'); version=version.slice(0,2); version=str2num(version);

//variables
md.qmu.variables.rheology_B=normal_uncertain('scaled_MaterialsRheologyB',1,.05);

//responses
md.qmu.responses.MaxVel=response_function('MaxVel',[],[0.0001 0.001 0.01 0.25 0.5 0.75 0.99 0.999 0.9999]);

//method
md.qmu.method     =dakota_method('nond_l');

//parameters
md.qmu.params.direct=true;
md.qmu.params.interval_type='forward';

if (version>=6,) {
	md.qmu.params.analysis_driver='matlab';
	md.qmu.params.evaluation_scheduling='master';
	md.qmu.params.processors_per_evaluation=2;
} else {
	md.qmu.params.analysis_driver='stressbalance';
	md.qmu.params.evaluation_concurrency=1;
}

//imperative! 
md.stressbalance.reltol=10^-10; //tighten for qmu analysese
md.qmu.isdakota=1;

//solve
md=solve(md,'Stressbalance','overwrite','y');

//Fields and tolerances to track changes
md.qmu.results=md.results.dakota;
md.results.dakota.importancefactors=importancefactors(md,'scaled_MaterialsRheologyB','MaxVel');
field_names     =['importancefactors'];
field_tolerances=[1e-10];
field_values=[
         md.results.dakota.importancefactors,
	];
