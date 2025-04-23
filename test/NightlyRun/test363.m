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
md.smb.lat = [1,2,3]';
md.smb.lon = [1,2,3]';
smbtime = [2,4,5];
Nlat = numel(md.smb.lat);
Nlon = numel(md.smb.lon);
Ntime = numel(smbtime);

md.smb.precipitation = [ones(Nlat*Nlon,Ntime)*0.1; smbtime];
md.smb.temperature = [reshape([1:Nlat*Nlon*Ntime], Nlat*Nlat, Ntime); smbtime];
%md.smb.temperature = [ones(Nlat*Nlon,Ntime)*3.6; smbtime];

md.transient.isthermal=0;
md.timestepping.cycle_forcing=0;
md=solve(md,'Transient');

