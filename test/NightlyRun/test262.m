%Test Name: SquareShelfTranIspddSicopolisMonthlyAnomaliesSSA2d
md=triangle(model(),'../Exp/Square.exp',150000.);
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareShelf.par');

% Use of SMBpddSicopolis
md.smb = SMBpddSicopolis();
% initalize pdd fields
md.smb=initialize(md.smb,md);
md.smb.s0p=md.geometry.surface;
md.smb.s0t=md.geometry.surface;

md.smb.monthlytemperatures=[]; md.smb.precipitation=[];
temp_ma_present=-10*ones(md.mesh.numberofvertices,1)-md.smb.rlaps*md.geometry.surface/1000;
temp_mj_present=10*ones(md.mesh.numberofvertices,1)-md.smb.rlaps*md.geometry.surface/1000;
precipitation=5*ones(md.mesh.numberofvertices,1);
for imonth=0:11
	md.smb.monthlytemperatures(1:md.mesh.numberofvertices,imonth+1)=md.materials.meltingpoint+temp_ma_present+(temp_mj_present-temp_ma_present)*sin(double(imonth+1-4)*pi/6.0);
	md.smb.precipitation(1:md.mesh.numberofvertices,imonth+1)=precipitation;
end

% time steps and resolution
md.timestepping.time_step=1;
md.settings.output_frequency=1;
md.timestepping.final_time=2;

md.transient.issmb=1;
md.transient.ismasstransport=1;
md.transient.isstressbalance=0;
md.transient.isthermal=0;

md.transient.requested_outputs={'default','TemperaturePDD'};
md.cluster=generic('name',oshostname(),'np',1); % 3 for the cluster

nv=md.mesh.numberofvertices;

%Run 1: legacy anomalies (one value per vertex, same for all 12 months)
md.smb.temperature_anomaly=1.5*ones(nv,1);
md.smb.precipitation_anomaly=1.1*ones(nv,1);
md=solve(md,'Transient');
res1=md.results.TransientSolution;

%Run 2: monthly anomalies (cell of 12 time series, time row in years last) holding the same constants.
%Must reproduce Run 1 exactly
md.smb.temperature_anomaly=cell(1,12);
md.smb.precipitation_anomaly=cell(1,12);
for imonth=1:12
	md.smb.temperature_anomaly{imonth}=[1.5*ones(nv,2); 0 2];
	md.smb.precipitation_anomaly{imonth}=[1.1*ones(nv,2); 0 2];
end
md=solve(md,'Transient');
res2=md.results.TransientSolution;

%Run 3: anomalies that vary with calendar month and with time
for imonth=1:12
	dT=2*sin(2*pi*(imonth-1)/12);
	md.smb.temperature_anomaly{imonth}=[dT*ones(nv,1) (dT+1)*ones(nv,1); 0 2];
	md.smb.precipitation_anomaly{imonth}=[(1+0.02*imonth)*ones(nv,1) (1.2+0.02*imonth)*ones(nv,1); 0 2];
end
md=solve(md,'Transient');
res3=md.results.TransientSolution;

%Consistency check: 12 identical constant monthly series must reproduce the legacy results
for i=1:2
	if max(abs(res1(i).SmbMassBalance-res2(i).SmbMassBalance))>1e-12*max(abs(res1(i).SmbMassBalance)) | ...
	   max(abs(res1(i).TemperaturePDD-res2(i).TemperaturePDD))>1e-12*max(abs(res1(i).TemperaturePDD))
		error('monthly anomalies do not reproduce the legacy (nv x 1) anomaly results');
	end
end

%Fields and tolerances to track changes
field_names     ={...
	'LegacyTemperaturePDD1','LegacySmbMassBalance1','LegacyTemperaturePDD2','LegacySmbMassBalance2',...
	'MonthlyConstTemperaturePDD1','MonthlyConstSmbMassBalance1','MonthlyConstTemperaturePDD2','MonthlyConstSmbMassBalance2',...
	'MonthlyVaryingTemperaturePDD1','MonthlyVaryingSmbMassBalance1','MonthlyVaryingTemperaturePDD2','MonthlyVaryingSmbMassBalance2'};
field_tolerances={1e-13,1e-13,1e-13,1e-13,...
	1e-13,1e-13,1e-13,1e-13,...
	1e-13,1e-13,1e-13,1e-13};
field_values={...
	(res1(1).TemperaturePDD),(res1(1).SmbMassBalance),(res1(2).TemperaturePDD),(res1(2).SmbMassBalance),...
	(res2(1).TemperaturePDD),(res2(1).SmbMassBalance),(res2(2).TemperaturePDD),(res2(2).SmbMassBalance),...
	(res3(1).TemperaturePDD),(res3(1).SmbMassBalance),(res3(2).TemperaturePDD),(res3(2).SmbMassBalance),...
	};
