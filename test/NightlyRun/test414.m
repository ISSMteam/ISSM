%Test Name: SquareSheetShelfDiadSSA3dDakotaMassFlux
md=triangle(model(),'../Exp/Square.exp',150000.);
md=setmask(md,'../Exp/SquareShelf.exp','');
md=parameterize(md,'../Par/SquareSheetShelf.par');
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',3);
md.materials.rho_ice=10^7; %involved in the mass flux, make it easy
md.geometry.thickness(:)=1; %make it easy
md.geometry.surface=md.geometry.base+md.geometry.thickness;

%constrain all velocities to 1 m/yr, in the y-direction
md.stressbalance.spcvx=zeros(md.mesh.numberofvertices,1);
md.stressbalance.spcvy=ones(md.mesh.numberofvertices,1);
md.stressbalance.spcvz=zeros(md.mesh.numberofvertices,1);

%Dakota options

%dakota version
version=IssmConfig('_DAKOTA_VERSION_'); version=version(1:3); version=str2num(version);

%partitioning
npart=20;
partition=partitioner(md,'package','chaco','npart',npart,'weighting','on')-1;

%variables
md.qmu.variables.drag_coefficient=normal_uncertain('descriptor','scaled_FrictionCoefficient',...
	'mean',ones(npart,1),...
	'stddev',.01*ones(npart,1),...
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

%method
md.qmu.method     =dakota_method('nond_l');

%parameters
md.qmu.params.direct=true;
md.qmu.params.interval_type='forward';
md.qmu.isdakota=1;
md.stressbalance.reltol=10^-5; %tighten for qmu analyses

if version>=6,
	md.qmu.params.analysis_driver='matlab';
	md.qmu.params.evaluation_scheduling='master';
	md.qmu.params.processors_per_evaluation=2;
else
	md.qmu.params.analysis_driver='stressbalance';
	md.qmu.params.evaluation_concurrency=1;
end


%solve
md=solve(md,'Stressbalance','overwrite','y');
md.qmu.results=md.results.dakota;

%Fields and tolerances to track changes
%ok, mass flux of 3 profiles should be -3 Gt/yr -3 Gt/yr and the sum, which is -6 Gt/yr
%we recover those mass fluxes through the mean of the response.
%also, we recover the max velo, which should be 1m/yr. 
%we put all that data in the moments, which we will use to test for success.
%also, check that the stddev are 0.
md.results.dakota.moments=[];
for i=1:8,
	md.results.dakota.moments=[md.results.dakota.moments md.results.dakota.dresp_out(i).mean];
end
for i=1:8,
	md.results.dakota.moments=[md.results.dakota.moments md.results.dakota.dresp_out(i).stddev];
end
field_names     ={'moments'};
field_tolerances={1e-11};
field_values={...
         md.results.dakota.moments,...
	};
