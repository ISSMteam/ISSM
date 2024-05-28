%Test Name: SquareShelfTranSemic
md=triangle(model(),'../Exp/Square.exp',150000.);
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareShelf.par');

% Use of SMBpddSicopolis
md.smb = SMBsemic();
% initalize pdd fields
md.smb=initialize(md.smb,md);
md.smb.s0gcm=md.geometry.surface;

ONES=ones(md.mesh.numberofvertices,1);
for iday=0:365
	md.smb.dailytemperature(1:md.mesh.numberofvertices,iday+1)=252.8739*ONES;
	md.smb.dailytemperature(md.mesh.numberofvertices+1,iday+1)=((iday+1)/12);
	md.smb.dailysnowfall(1:md.mesh.numberofvertices,iday+1)=8.5503e-09*ONES;
	md.smb.dailysnowfall(md.mesh.numberofvertices+1,iday+1)=((iday+1)/12);
	md.smb.dailyrainfall(1:md.mesh.numberofvertices,iday+1)=1.7296e-09*ONES;
	md.smb.dailyrainfall(md.mesh.numberofvertices+1,iday+1)=((iday+1)/12);
	md.smb.dailydsradiation(1:md.mesh.numberofvertices,iday+1)=128.1702*ONES;
	md.smb.dailydsradiation(md.mesh.numberofvertices+1,iday+1)=((iday+1)/12);
	md.smb.dailydlradiation(1:md.mesh.numberofvertices,iday+1)=176.5667*ONES;
	md.smb.dailydlradiation(md.mesh.numberofvertices+1,iday+1)=((iday+1)/12);
	md.smb.dailywindspeed(1:md.mesh.numberofvertices,iday+1)=6.0741*ONES;
	md.smb.dailywindspeed(md.mesh.numberofvertices+1,iday+1)=((iday+1)/12);
	md.smb.dailyairdensity(1:md.mesh.numberofvertices,iday+1)=1.0729*ONES;
	md.smb.dailyairdensity(md.mesh.numberofvertices+1,iday+1)=((iday+1)/12);
	md.smb.dailyairhumidity(1:md.mesh.numberofvertices,iday+1)=9.6667e-04*ONES;
	md.smb.dailyairhumidity(md.mesh.numberofvertices+1,iday+1)=((iday+1)/12); 
	md.smb.dailypressure(1:md.mesh.numberofvertices,iday+1)=7.7841e+04*ONES;
	md.smb.dailypressure(md.mesh.numberofvertices+1,iday+1)=((iday+1)/12); 
end

% time steps and resolution
md.timestepping.time_step=0.5;
md.settings.output_frequency=1;
md.timestepping.final_time=1;

md.transient.issmb=1;
md.transient.ismasstransport=0;
md.transient.isstressbalance=0;
md.transient.isthermal=0;

md.transient.requested_outputs={'default','TemperatureSEMIC'};
md.cluster=generic('name',oshostname(),'np',4); % 3 for the cluster
md=solve(md,'Transient');

%Fields and tolerances to track changes
field_names     ={'TemperatureSEMIC1','SmbMassBalance1','TemperatureSEMIC2','SmbMassBalance2'};
field_tolerances={1e-13,1e-13,1e-13,1e-13};
field_values={...
	(md.results.TransientSolution(1).TemperatureSEMIC),...
	(md.results.TransientSolution(1).SmbMassBalance),...
	(md.results.TransientSolution(2).TemperatureSEMIC),...
	(md.results.TransientSolution(2).SmbMassBalance),...
	};
