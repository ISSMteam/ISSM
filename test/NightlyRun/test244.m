%Test Name: SquareShelfSMBGembDakota
md=triangle(model(),'../Exp/Square.exp',300000.);
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareShelf.par');
md=setflowequation(md,'SSA','all');
md.materials.rho_ice=910;
md.cluster=generic('name',oshostname(),'np',3);
md.geometry.bed=md.geometry.base;

% Use of Gemb method for SMB computation
md.smb = SMBgemb(md.mesh);
md.smb.dsnowIdx = 0;
md.smb.swIdx = 1;

%load hourly surface forcing date from 1979 to 2009:
inputs=load('../Data/gemb_input.mat');

%setup the inputs:
md.smb.Ta=[repmat(inputs.Ta0',md.mesh.numberofelements,1);inputs.dateN'];
md.smb.V=[repmat(inputs.V0',md.mesh.numberofelements,1);inputs.dateN'];
md.smb.dswrf=[repmat(inputs.dsw0',md.mesh.numberofelements,1);inputs.dateN'];
md.smb.dlwrf=[repmat(inputs.dlw0',md.mesh.numberofelements,1);inputs.dateN'];
md.smb.P=[repmat(inputs.P0',md.mesh.numberofelements,1);inputs.dateN'];
md.smb.eAir=[repmat(inputs.eAir0',md.mesh.numberofelements,1);inputs.dateN'];
md.smb.pAir=[repmat(inputs.pAir0',md.mesh.numberofelements,1);inputs.dateN'];
md.smb.Vz=repmat(inputs.LP.Vz,md.mesh.numberofelements,1);
md.smb.Tz=repmat(inputs.LP.Tz,md.mesh.numberofelements,1);
md.smb.Tmean=repmat(inputs.LP.Tmean,md.mesh.numberofelements,1);
md.smb.C=repmat(inputs.LP.C,md.mesh.numberofelements,1);

%smb settings
md.smb.requested_outputs={'SmbDz','SmbT','SmbD','SmbRe','SmbGdn','SmbGsp','SmbEC','SmbA','SmbMassBalance'};

%only run smb core:
md.transient.isstressbalance=0;
md.transient.ismasstransport=1;
md.transient.isthermal=0;

%time stepping:
md.timestepping.start_time=1965;
md.timestepping.final_time=1965.75;
md.timestepping.time_step=1/365.0;
md.timestepping.interp_forcing=0;

%dakota version
version=IssmConfig('_DAKOTA_VERSION_'); version=version(1:3); version=str2num(version);

%partitioning
partition=partitioner(md,'package','linear','type','element','npart',md.mesh.numberofelements)-1;

%variables
md.qmu.variables.surface_mass_balanceC=normal_uncertain('descriptor','scaled_SmbC',...
	'mean',1*ones(md.mesh.numberofelements,1),...
	'stddev',.5*ones(md.mesh.numberofelements,1),...
	'partition',partition);

Tmin=273;
telms=min(md.smb.Ta(1:end-1,:),[],2);
mint_on_partition=telms;
for pa=1:length(mint_on_partition)
	vi=find(partition+1 == pa);
	mint=telms(vi).*1.05;
	pos=find(mint < Tmin);
	mint(pos)=Tmin;
	mint_on_partition(pa)=max(mint./telms(vi));
end
mint_on_partition(isnan(mint_on_partition)) = 10^-10;
md.qmu.variables.surface_mass_balanceTa=uniform_uncertain('descriptor','scaled_SmbTa',...
	'lower',.95*ones(md.mesh.numberofelements,1),...
	'upper',max(min(max(1.05,mint_on_partition),0.9999),0.0001),...
	'partition',partition);

%responses
md.qmu.responses.IceVolume=response_function('descriptor','IceVolume');
md.qmu.responses.IceMass=response_function('descriptor','IceMass');
md.qmu.responses.TotalSmb=response_function('descriptor','TotalSmb');

%  nond_sampling study
md.qmu.method=dakota_method('nond_samp');
md.qmu.method(end)=dmeth_params_set(md.qmu.method(end),'seed',1234,'samples',3,'sample_type','lhs');
dver=textscan(IssmConfig('_DAKOTA_VERSION_'),'%[0123456789].%[0123456789].%[0123456789]');
if ((str2num(dver{1}{1})==4 && str2num(dver{2}{1})>2) || str2num(dver{1}{1})>4)
	md.qmu.method(end)=dmeth_params_set(md.qmu.method(end),'rng','rnum2');
end

%parameters
md.qmu.params.direct=true;
md.qmu.params.analysis_components='';
md.qmu.params.interval_type='forward';
md.qmu.params.tabular_graphics_data=true;
md.qmu.isdakota=1;

if version>=6,
	md.qmu.params.analysis_driver='matlab';
	md.qmu.params.evaluation_scheduling='master';
	md.qmu.params.processors_per_evaluation=2;
else
	md.qmu.params.analysis_driver='stressbalance';
	md.qmu.params.evaluation_concurrency=1;
end

md.stressbalance.reltol=10^-5; %tighten for qmu analyses
md.transient.requested_outputs={'IceVolume','TotalSmb','IceMass'};

%solve
md=solve(md,'Transient','overwrite','y');
md.qmu.results=md.results.dakota;

%Fields and tolerances to track changes
md.results.dakota.moments=[];
for i=1:3,
	md.results.dakota.moments=[md.results.dakota.moments md.results.dakota.dresp_out(i).mean];
end
for i=1:3,
	md.results.dakota.moments=[md.results.dakota.moments md.results.dakota.dresp_out(i).stddev];
end
field_names     ={'moments'};
field_tolerances={2e-6};
field_values={...
	md.results.dakota.moments,...
	};

