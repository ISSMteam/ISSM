%Test Name: SquareShelfTranIspddSicopolisSSA2d
md=triangle(model(),'../Exp/Square.exp',150000.);
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareShelf.par');

% Use of SMBpddSicopolis
md.smb = SMBpddSicopolis();
% initalize pdd fields
md.smb=initialize(md.smb,md);
md.smb.s0p=md.geometry.surface;
md.smb.s0t=md.geometry.surface;

% 
md.smb.monthlytemperatures=[]; md.smb.precipitation=[]; md.smb.precipitation=[];
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
md=solve(md,'Transient');

%Fields and tolerances to track changes
field_names     ={'TemperaturePDD1','SmbMassBalance1','TemperaturePDD2','SmbMassBalance2'};
field_tolerances={1e-13,1e-13,1e-13,1e-13};
field_values={...
	(md.results.TransientSolution(1).TemperaturePDD),...
	(md.results.TransientSolution(1).SmbMassBalance),...
	(md.results.TransientSolution(2).TemperaturePDD),...
	(md.results.TransientSolution(2).SmbMassBalance),...
	};
