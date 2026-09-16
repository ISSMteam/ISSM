%Test Name: SquareShelfTranSemicStandaloneTransient

% Initialize coordinates from '../Data/semic_transect_input.txt'
md=model;
md.miscellaneous.name='SquareShelfTranSemicStandaloneTransient';
xy = [[0,1,0,1,0,1,0];...
		[0,0,1,1,2,2,3]];
elements= [[2,3,1];...
	[2,4,3];...
	[3,4,5];...
	[4,5,6];...
	[5,6,7]];
segments=[[1,2,1];...
	[1,3,1];...
	[2,4,2];...
	[3,5,3];...
	[4,6,4];...
	[5,7,5];...
	[6,7,5]];
md.mesh.x=xy(1,:)';
md.mesh.y=xy(2,:)';
md.mesh.elements=elements;
md.mesh.segments=segments;
md.mesh.numberofvertices=length(md.mesh.x);
md.mesh.numberofelements=length(elements);
nV=md.mesh.numberofvertices;

% Set initial settings
md.geometry.surface=10*ones(nV,1);
md.geometry.base   =zeros(nV,1);
md.geometry.bed    =zeros(nV,1);
md.geometry.thickness=md.geometry.surface-md.geometry.base;

% Set mask
md.mask.ice_levelset=-1*ones(nV,1);
md.mask.ocean_levelset=ones(nV,1);

% Use of SMBsemic
md.smb = SMBsemic();

% Initialize parameters
md.smb.ismethod=1; % Use real transient mode in SEMIC
md.smb.albedo_scheme=2; % Denby albedo scheme

md.smb.alb_smax=0.81;
md.smb.alb_smin=0.6;
md.smb.tmax=273.15;
md.smb.tmin=273.15-10;
md.smb.hcrit=0.028;
md.smb.rcrit=0.79;
md.smb.Tamp=3.5*ones(nV,1);
md.smb.mask=[1,1,1,2,2,2,2]';
md.smb.s0gcm      = md.geometry.surface;

% Initialize SEMIC values
md.initialization.temperature=(273.15-10)*ones(nV,1);
md.smb.albedo     = 0.8*ones(nV,1);
md.smb.albedo_snow= 0.8*ones(nV,1);
md.smb.hsnow      = 5*ones(nV,1);
md.smb.hice       = zeros(nV,1); 
md.smb.qmr        = zeros(nV,1);

% Load forcing dataset
% Number of stations : 7
% Number of days: 365
nx=7;
times=[1:365];
temp=readtable('../Data/semic_transect_input.txt');
forc=struct;
forc.dailysnowfall   =transpose(table2array(temp(:,0*nx+1:1*nx)));
forc.dailyrainfall   =transpose(table2array(temp(:,1*nx+1:2*nx)));
forc.dailydsradiation=transpose(table2array(temp(:,2*nx+1:3*nx)));
forc.dailydlradiation=transpose(table2array(temp(:,3*nx+1:4*nx)));
forc.dailywindspeed  =transpose(table2array(temp(:,4*nx+1:5*nx)));
forc.dailypressure   =transpose(table2array(temp(:,5*nx+1:6*nx)));
forc.dailyairdensity =transpose(table2array(temp(:,6*nx+1:7*nx)));
forc.dailyairhumidity=transpose(table2array(temp(:,7*nx+1:8*nx)));
forc.dailytemperature=transpose(table2array(temp(:,8*nx+1:9*nx)));

fields=fieldnames(forc);
for i=1:length(fields)
	field=fields{i};
	md.smb.(field) = [forc.(field); times];
end

% time steps and resolution
md.timestepping.time_step=1/365; % 1 day
md.settings.output_frequency=1;
md.timestepping.final_time=1;

md.transient.issmb=1;
md.transient.ismasstransport=0;
md.transient.isstressbalance=0;
md.transient.isthermal=0;

md.transient.requested_outputs={'default','TemperatureSEMIC','SmbMelt'};
md.cluster=generic('name',oshostname(),'np',1); % 3 for the cluster
md.verbose.solution=0;
md=solve(md,'Transient');

%Fields and tolerances to track changes
field_names     ={'TemperatureSEMIC1','SmbMassBalance1','SmbMelt1','TemperatureSEMIC2','SmbMassBalance2','SmbMelt2'};
field_tolerances={1e-13,1e-13,1e-13,1e-13,1e-13,1e-13};
field_values={...
	(md.results.TransientSolution(1).TemperatureSEMIC),...
	(md.results.TransientSolution(1).SmbMassBalance),...
	(md.results.TransientSolution(1).SmbMelt),...
	(md.results.TransientSolution(365).TemperatureSEMIC),...
	(md.results.TransientSolution(365).SmbMassBalance),...
	(md.results.TransientSolution(365).SmbMelt),...
	};
