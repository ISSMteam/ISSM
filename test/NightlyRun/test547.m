%Test Name: PigTranFrontalforcingsrignotRunoff
md=triangle(model(),'../Exp/Pig.exp',10000.);
md=setmask(md,'../Exp/PigShelves.exp','../Exp/PigIslands.exp');
md=parameterize(md,'../Par/Pig.par');
md=setflowequation(md,'SSA','all');
md.timestepping.time_step=1;
md.timestepping.final_time=3;

%separate domain in 2 basins
idbasin = zeros(md.mesh.numberofelements,1);
iid1    = find(md.mesh.x<=-1.6e6);
for ii=1:md.mesh.numberofelements
	for vertex=1:3
		if any(iid1==md.mesh.elements(ii,vertex)) %one vertex in basin 1
			idbasin(ii) = 1;
		end
	end
	if idbasin(ii)==0 %no vertex was found in basin 1
		idbasin(ii) = 2;
	end
end

%SMB: PDD scheme providing the runoff that drives the subglacial discharge
%(monthly cycle of -18C to +2C at sea level, 0.4 m/yr precip: ~0.7-1 m/yr of runoff on the low-lying shelves, none at high elevation)
nv=md.mesh.numberofvertices;
md.smb = SMBpddSicopolis();
md.smb=initialize(md.smb,md);
surface=max(md.geometry.surface,0);
md.smb.s0p=surface;
md.smb.s0t=surface;
temp_ma_present=-18*ones(nv,1)-md.smb.rlaps*surface/1000;
temp_mj_present= 2*ones(nv,1)-md.smb.rlaps*surface/1000;
precipitation=0.4*ones(nv,1);
md.smb.monthlytemperatures=zeros(nv,12);
md.smb.precipitation=zeros(nv,12);
for imonth=0:11
	md.smb.monthlytemperatures(1:nv,imonth+1)=md.materials.meltingpoint+temp_ma_present+(temp_mj_present-temp_ma_present)*sin(double(imonth+1-4)*pi/6.0);
	md.smb.precipitation(1:nv,imonth+1)=precipitation;
end

%calving and frontalforcings parameters
md.mask.ice_levelset = 1e4*(md.mask.ice_levelset + 0.5);
md.calving=calvingvonmises();
md.frontalforcings=frontalforcingsrignot();
md.frontalforcings.num_basins = 2;
md.frontalforcings.basin_id = idbasin;
md.frontalforcings.isdischargerunoff = 1; %discharge = basin-integrated SMB runoff, subglacial_discharge is not used
md.frontalforcings.thermalforcing = 0.5*ones(nv,1);
for elem=1:md.mesh.numberofelements
	if idbasin(elem)==2
		md.frontalforcings.thermalforcing(md.mesh.elements(elem,1:3)) = 1.5;
	end
end

md.transient.issmb = 1;
md.transient.ismovingfront = 1;
md.levelset.spclevelset = NaN(nv,1);
md.levelset.migration_max = 1e10;

md.transient.requested_outputs = {'default','CalvingMeltingrate','SmbRunoff'};
md.cluster=generic('name',oshostname(),'np',2);
md=solve(md,'Transient');

%Fields and tolerances to track changes
field_names ={...
	'Vel1','Thickness1','MaskIceLevelset1','CalvingMeltingrate1','SmbRunoff1',...
	'Vel2','Thickness2','MaskIceLevelset2','CalvingMeltingrate2','SmbRunoff2',...
	'Vel3','Thickness3','MaskIceLevelset3','CalvingMeltingrate3','SmbRunoff3',...
	};
field_tolerances={...
	2e-11,9e-12,1e-11,1e-9,1e-13,...
	2e-11,2e-11,1e-11,1e-9,2e-13,...
	1e-10,1e-10,1e-9,1e-9,2e-13,...
	};
field_values={...
	(md.results.TransientSolution(1).Vel),...
	(md.results.TransientSolution(1).Thickness),...
	(md.results.TransientSolution(1).MaskIceLevelset),...
	(md.results.TransientSolution(1).CalvingMeltingrate),...
	(md.results.TransientSolution(1).SmbRunoff),...
	(md.results.TransientSolution(2).Vel),...
	(md.results.TransientSolution(2).Thickness),...
	(md.results.TransientSolution(2).MaskIceLevelset),...
	(md.results.TransientSolution(2).CalvingMeltingrate),...
	(md.results.TransientSolution(2).SmbRunoff),...
	(md.results.TransientSolution(3).Vel),...
	(md.results.TransientSolution(3).Thickness),...
	(md.results.TransientSolution(3).MaskIceLevelset),...
	(md.results.TransientSolution(3).CalvingMeltingrate),...
	(md.results.TransientSolution(3).SmbRunoff),...
	};
