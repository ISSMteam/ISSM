//Test Name: SquareShelfTranForceNeg2dDakotaSamp
var md = new model();
triangle(md,square[0],180000.);
setmask(md,'all','');
parameterize(md);
setflowequation(md,'SSA','all');
//md.cluster=generic('name',oshostname(),'np',3);

md.timestepping.time_step=1;
md.settings.output_frequency=1;
md.timestepping.final_time=4;

smb=ones(md.mesh.numberofvertices,1);
for (var i = 0; i < smb.length; ++i) {
    smb[i][0] *= 3.6;
    smb[i].push(smb[i][0]*-1);
}

md.smb.mass_balance= smb.slice();
md.smb.mass_balance[md.smb.mass_balance.length-1] = [1.5, 3.];
md.trans.isthermal=0;
//Dakota options

//dakota version
version=IssmConfig('_DAKOTA_VERSION_'); version=version.toString().slice(0,2);

//partitioning
npart=20;
partitioner(md,'package','chaco','npart',npart,'weighting','on');
md.qmu.partition=md.qmu.partition-1;

//variables
md.qmu.variables.surface_mass_balance=normal_uncertain('scaled_SmbMassBalance',1,0.1);

//responses
md.qmu.responses.MaxVel=response_function('MaxVel',[],[0.0001, 0.001, 0.01, 0.25, 0.5, 0.75, 0.99, 0.999, 0.9999]);
md.qmu.responses.IceVolume=response_function('IceVolume',[],[0.0001, 0.001, 0.01, 0.25, 0.5, 0.75, 0.99, 0.999, 0.9999]);
md.qmu.responses.MassFlux1=response_function('indexed_MassFlux_1',[],[0.0001, 0.001, 0.01, 0.25, 0.5, 0.75, 0.99, 0.999, 0.9999]);
md.qmu.responses.MassFlux2=response_function('indexed_MassFlux_2',[],[0.0001, 0.001, 0.01, 0.25, 0.5, 0.75, 0.99, 0.999, 0.9999]);
md.qmu.responses.MassFlux3=response_function('indexed_MassFlux_3',[],[0.0001, 0.001, 0.01, 0.25, 0.5, 0.75, 0.99, 0.999, 0.9999]);
md.qmu.responses.MassFlux4=response_function('indexed_MassFlux_4',[],[0.0001, 0.001, 0.01, 0.25, 0.5, 0.75, 0.99, 0.999, 0.9999]);
md.qmu.responses.MassFlux5=response_function('indexed_MassFlux_5',[],[0.0001, 0.001, 0.01, 0.25, 0.5, 0.75, 0.99, 0.999, 0.9999]);
md.qmu.responses.massFlux6=response_function('indexed_MassFlux_6',[],[0.0001, 0.001, 0.01, 0.25, 0.5, 0.75, 0.99, 0.999, 0.9999]);

//mass flux profiles
md.qmu.mass_flux_profiles=['../Exp/MassFlux1.exp','../Exp/MassFlux2.exp','../Exp/MassFlux3.exp','../Exp/MassFlux4.exp','../Exp/MassFlux5.exp','../Exp/MassFlux6.exp'];
md.qmu.mass_flux_profile_directory=pwd;

////  nond_sampling study
md.qmu.method=dakota_method('nond_samp');
md.qmu.method[md.qmu.method.length-1]=dmeth_params_set(md.qmu.method[md.qmu.method.length-1],'seed',1234,'samples',20,'sample_type','lhs');
dver=textscan(IssmConfig('_DAKOTA_VERSION_'),'//[0123456789].//[0123456789].//[0123456789]');
if (((str2num(dver[1][1])==4 && str2num(dver[2][1])>2) || str2num(dver[1][1])>4)) {
	md.qmu.method[md.qmu.method.length-1]=dmeth_params_set(md.qmu.method(end),'rng','rnum2');
}

//parameters
md.qmu.params.direct=true;
md.qmu.params.analysis_components='';
md.qmu.params.interval_type='forward';
md.qmu.params.tabular_graphics_data=true;
md.qmu.isdakota=1;

if (version>=6) {
	md.qmu.params.analysis_driver='matlab';
	md.qmu.params.evaluation_scheduling='master';
	md.qmu.params.processors_per_evaluation=2;
} else {
	md.qmu.params.analysis_driver='stressbalance';
	md.qmu.params.evaluation_concurrency=1;
}


md.stressbalance.reltol=Math.pow(10,-5); //tighten for qmu analyses
md.trans.requested_outputs=['IceVolume'];

//solve
md=solve(md,'Transient','overwrite','y');
md.qmu.results=md.results.dakota;

//Fields and tolerances to track changes
md.results.dakota.moments=[];
for (var i = 0; i < 8; ++i) {
	md.results.dakota.moments.push(md.results.dakota.dresp_out[i].mean);
}
for (var i = 0; i < 8; ++i) {
	md.results.dakota.moments.push(md.results.dakota.dresp_out[i].stddev);
}
field_names     =['moments'];
field_tolerances=[1e-11];
field_values=[
         md.results.dakota.moments,
	];
