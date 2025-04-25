%Test Name: SquareSheetConstrainedSmbpddGCM2d
md=triangle(model(),'../Exp/Square.exp',150000.);
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareShelf.par');
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',3);

md.timestepping.time_step=0.5;
md.settings.output_frequency=1;
md.timestepping.final_time=6.;

%Set up SMB
md.smb=SMBpddGCM();
md.smb = md.smb.initialize(md);
md.smb.x_grid = linspace(0,1e6,3)';
md.smb.y_grid = linspace(0,1e6,4)';
smbtime = [2,4,5];
Nx = numel(md.smb.x_grid);
Ny = numel(md.smb.y_grid);
Ntime = numel(smbtime);

md.smb.precipitation = [ones(Nx*Ny,Ntime)*0.1; smbtime];
md.smb.temperature = [reshape([1:Nx*Ny*Ntime], Nx*Ny, Ntime)+273.15-(Nx*Ny);smbtime];
%md.smb.temperature = [ones(Nlat*Nlon,Ntime)*3.6; smbtime];

md.transient.isthermal=0;
md.timestepping.cycle_forcing=0;

md.transient.requested_outputs={'default','SmbTemperature', 'SmbPrecipitation','SmbAccumulation','SmbAblation','SmbMeanTemperature'};
%md.verbose = verbose('solution',1);
md=solve(md,'Transient');


t = [md.results.TransientSolution(:).SmbTemperature];
p = [md.results.TransientSolution(:).SmbPrecipitation];
meant = [md.results.TransientSolution(:).SmbMeanTemperature];
