%Test Name: ThermalEnthBasalcondsTrans
md=triangle(model(),'../Exp/Square.exp',300000.);
md=setmask(md,'','');
md=parameterize(md,'../Par/SquareThermal.par');

h=100.;
md.geometry.thickness=h*ones(md.mesh.numberofvertices,1);
md.geometry.base=-h*ones(md.mesh.numberofvertices,1);
md.geometry.surface=md.geometry.base+md.geometry.thickness;

md=extrude(md,41,2.);
md=setflowequation(md,'HO','all');
md.thermal.isenthalpy=1;
md.thermal.isdynamicbasalspc=1;

%Basal forcing
Ts=273.15-3.; Tb=273.15-1.; Tsw=Tb;
qgeo=md.materials.thermalconductivity/max(md.geometry.thickness)*(Tb-Ts); %qgeo=kappa*(Tb-Ts)/H
md.basalforcings.geothermalflux(find(md.mesh.vertexonbase))=qgeo;
md.initialization.temperature=qgeo/md.materials.thermalconductivity.*(md.geometry.surface-md.mesh.z)+Ts;

%Surface forcing
pos=find(md.mesh.vertexonsurface);
SPC_cold=NaN(md.mesh.numberofvertices,1);
SPC_warm=NaN(md.mesh.numberofvertices,1);
SPC_cold(pos)=Ts;
SPC_warm(pos)=Tsw;
md.thermal.spctemperature=SPC_cold;
md.timestepping.time_step=5.;
t0=0.;
t1=10.;
t2=100.;
md.timestepping.final_time=400.;
md.thermal.spctemperature=[SPC_cold SPC_cold SPC_warm SPC_warm SPC_cold];
md.thermal.spctemperature=[md.thermal.spctemperature; t0 t1-1 t1 t2-1 t2];

%Additional settings
md.transient.ismasstransport=0;
md.transient.isstressbalance=0;
md.transient.issmb=1;
md.transient.isthermal=1;
md.thermal.stabilization = 0;

%Go solve
md.verbose=verbose(0);
md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Transient');

%Fields and tolerances to track changes
field_names     ={'Enthalpy1','Temperature1','Waterfraction1','BasalMeltingRate1','Watercolumn1',...
						'Enthalpy2','Temperature2','Waterfraction2','BasalMeltingRate2','Watercolumn2',...
						'Enthalpy3','Temperature3','Waterfraction3','BasalMeltingRate3','Watercolumn3',...
						'Enthalpy4','Temperature4','Waterfraction4','BasalMeltingRate4','Watercolumn4'};
field_tolerances={1.e-10,1.e-10,1.e-10,1.e-9,1.e-10,...
						1.e-10,1.e-10,1.e-10,2.e-9,2.e-10,...
						1.e-10,1.e-10,1.e-10,2.e-9,1.e-10,...
						1.e-10,1.e-10,1.e-10,2.e-9,1.e-10};
i1=1;	i2=ceil(t2/md.timestepping.time_step)+2;	i3=ceil(md.timestepping.final_time/(2.*md.timestepping.time_step));	i4=size(md.results.TransientSolution,2);
field_values={...
	(md.results.TransientSolution(i1).Enthalpy),...
	(md.results.TransientSolution(i1).Temperature),...
	(md.results.TransientSolution(i1).Waterfraction),...
	(md.results.TransientSolution(i1).BasalforcingsGroundediceMeltingRate),...
	(md.results.TransientSolution(i1).Watercolumn),...
	(md.results.TransientSolution(i2).Enthalpy),...
	(md.results.TransientSolution(i2).Temperature),...
	(md.results.TransientSolution(i2).Waterfraction),...
	(md.results.TransientSolution(i2).BasalforcingsGroundediceMeltingRate),...
	(md.results.TransientSolution(i2).Watercolumn),...
	(md.results.TransientSolution(i3).Enthalpy),...
	(md.results.TransientSolution(i3).Temperature),...
	(md.results.TransientSolution(i3).Waterfraction),...
	(md.results.TransientSolution(i3).BasalforcingsGroundediceMeltingRate),...
	(md.results.TransientSolution(i3).Watercolumn),...
	(md.results.TransientSolution(i4).Enthalpy),...
	(md.results.TransientSolution(i4).Temperature),...
	(md.results.TransientSolution(i4).Waterfraction),...
	(md.results.TransientSolution(i4).BasalforcingsGroundediceMeltingRate),...
	(md.results.TransientSolution(i4).Watercolumn),...
	};
