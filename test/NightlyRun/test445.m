%Test Name: SquareSheetShelfSteaEnthalpyHO3dDakotaSampNeff

% TODO:
% - Figure out why test fails intermittently on Mac with "Index exceeds array 
% bounds."
%

md=triangle(model(),'../Exp/Square.exp',150000.);
md=setmask(md,'../Exp/SquareShelf.exp','');
md=parameterize(md,'../Par/SquareSheetShelf.par');
md=extrude(md,3,2.);
md=setflowequation(md,'HO','all');
md.cluster=generic('name',oshostname(),'np',3);
md.timestepping.time_step=0.;
md.thermal.isenthalpy=1;
md.thermal.isdynamicbasalspc=1;
md.initialization.waterfraction=zeros(md.mesh.numberofvertices,1);
md.initialization.watercolumn=zeros(md.mesh.numberofvertices,1);

md.friction.coupling=3;
md.friction.effective_pressure=md.materials.rho_ice*md.constants.g*md.geometry.thickness+md.materials.rho_water*md.constants.g*md.geometry.base;

%dakota version
version=IssmConfig('_DAKOTA_VERSION_');
version=version(1:3);
version=str2num(version);

%partitioning
npart=10;
partition=partitioner(md,'package','chaco','npart',npart,'weighting','on')-1;
md.qmu.isdakota=1;

%variables
md.qmu.variables.neff=normal_uncertain('descriptor','scaled_FrictionEffectivePressure',...
	'mean',ones(npart,1),...
	'stddev',.05*ones(npart,1),...
	'partition',partition);
md.qmu.variables.geoflux=normal_uncertain('descriptor','scaled_BasalforcingsGeothermalflux',...
	'mean',ones(npart,1),...
	'stddev',.05*ones(npart,1),...
	'partition',partition);

%responses
md.qmu.responses.MaxVel=response_function('descriptor','MaxVel');
md.qmu.responses.MassFlux1=response_function('descriptor','indexed_MassFlux_1');
md.qmu.responses.MassFlux2=response_function('descriptor','indexed_MassFlux_2');
md.qmu.responses.MassFlux3=response_function('descriptor','indexed_MassFlux_3');
md.qmu.responses.MassFlux4=response_function('descriptor','indexed_MassFlux_4');
md.qmu.responses.MassFlux5=response_function('descriptor','indexed_MassFlux_5');
md.qmu.responses.MassFlux6=response_function('descriptor','indexed_MassFlux_6');
md.qmu.responses.MassFlux7=response_function('descriptor','indexed_MassFlux_7');

%mass flux profiles
md.qmu.mass_flux_profiles={'../Exp/MassFlux1.exp','../Exp/MassFlux2.exp','../Exp/MassFlux3.exp','../Exp/MassFlux4.exp','../Exp/MassFlux5.exp','../Exp/MassFlux6.exp','../Exp/Square.exp'};
md.qmu.mass_flux_profile_directory=pwd;

md.qmu.method     =dakota_method('nond_samp');
md.qmu.method(end)=dmeth_params_set(md.qmu.method(end),...
	'seed',1234,...
	'samples',20,...
	'sample_type','random');

%%  a variety of parameters
md.qmu.params.direct=true;
md.qmu.params.analysis_components='';
md.qmu.params.tabular_graphics_data=true;

if version>=6,
	md.qmu.params.analysis_driver='matlab';
	md.qmu.params.evaluation_scheduling='master';
	md.qmu.params.processors_per_evaluation=2;
else
	md.qmu.params.analysis_driver='stressbalance';
	md.qmu.params.evaluation_concurrency=1;
end


md.stressbalance.reltol=10^-5; %tighten for qmu analyses

md=solve(md,'Steadystate','overwrite','y');

%Fields and tolerances to track changes
md.qmu.results=md.results.dakota;

%we put all the mean and stdev data in the montecarlo field, which we will use to test for success.
md.results.dakota.montecarlo=[];
for i=1:8,
	md.results.dakota.montecarlo=[md.results.dakota.montecarlo md.results.dakota.dresp_out(i).mean];
end
for i=1:8,
	md.results.dakota.montecarlo=[md.results.dakota.montecarlo md.results.dakota.dresp_out(i).stddev];
end
field_names     ={'montecarlo'};
field_tolerances={2e-10};
field_values={...
	md.results.dakota.montecarlo,...
	};

