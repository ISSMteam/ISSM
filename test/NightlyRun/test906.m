%Test Name: SquareNoDynHydrologyDCtimeadapt
md=triangle(model(), '../Exp/Square.exp', 200000.);
md=setmask(md, '', '');
md=parameterize(md, '../Par/SquareNoDyn.par');
md.cluster=generic('name', oshostname(), 'np', 1);

md.transient=deactivateall(md.transient);
md.transient.ishydrology=1;
md.hydrology=(hydrologydc);
md.hydrology=initialize(md.hydrology,md);

%Hydro Model Parameters
md.hydrology.isefficientlayer=1;
md.hydrology.sedimentlimit_flag=1;
md.hydrology.sedimentlimit=500;
md.hydrology.rel_tol=1.0e-6;
md.hydrology.penalty_lock=10;
md.hydrology.max_iter=100;
md.hydrology.transfer_flag=1;
md.hydrology.unconfined_flag=0;
md.hydrology.leakage_factor=5.0e-10;
md.hydrology.mask_thawed_node=ones(md.mesh.numberofvertices,1);

%Sediment Parameters
md.hydrology.sediment_porosity=0.4;
md.hydrology.sediment_transmitivity=3.0*ones(md.mesh.numberofvertices,1);
md.hydrology.sediment_thickness=10.;

%Epl Parameters
md.hydrology.mask_eplactive_node=zeros(md.mesh.numberofvertices,1);
md.hydrology.epl_conductivity=30.;
md.hydrology.epl_initial_thickness=1.0;
md.hydrology.epl_colapse_thickness=1.0e-3;
md.hydrology.epl_thick_comp=0;
md.hydrology.epl_max_thickness=10;
md.hydrology.eplflip_lock=0;

%Initialisation
md.initialization.sediment_head=md.geometry.base;
md.initialization.epl_head=md.geometry.base;
md.initialization.epl_thickness=ones(md.mesh.numberofvertices,1);

%Boundary conditions
md.hydrology.spcsediment_head=NaN*ones(md.mesh.numberofvertices,1);
md.hydrology.spcepl_head=NaN*ones(md.mesh.numberofvertices,1);

%Forcing
stepping=3600./md.constants.yts;
endtime=0.3;
times=[0.:stepping:endtime];
input=find(abs(times-0.3) > 0.2);
inputval=times*0+100;
inputval(input)=0;

%basal forcing inputval
md.basalforcings.groundedice_melting_rate=repmat(inputval,md.mesh.numberofvertices,1);
md.basalforcings.groundedice_melting_rate=[[md.basalforcings.groundedice_melting_rate];[times]];

%Time
md.timestepping.final_time=endtime;
md.settings.output_frequency=1.0;

md.timestepping.average_forcing=1;
md.hydrology.step_adapt=1;
md.hydrology.steps_per_step=1;
md.timestepping.time_step=0.1;


md=solve(md, 'Transient');

field_names={'SedimentWaterHead2', 'EplWaterHead2', 'SedimentHeadResidual2',...
	     'SedimentHeadSubstep2', 'EplHeadSubstep2', 'HydrologySubsteps2', 'HydrologySubTime2',...
	     'SedimentWaterHead3', 'EplWaterHead3', 'SedimentHeadResidual3'};
field_tolerances={1e-13, 1e-13, 1e-13,...
		  1e-13, 1e-13, 1e-13, 1e-13,...
		  1e-13, 1e-13, 1e-13};
field_values={md.results.TransientSolution(2).SedimentHead,...
	      md.results.TransientSolution(2).EplHead,...
	      md.results.TransientSolution(2).SedimentHeadResidual,...
	      md.results.TransientSolution(2).SedimentHeadSubstep,...
	      md.results.TransientSolution(2).EplHeadSubstep,...
	      md.results.TransientSolution(2).HydrologySubsteps,...
	      md.results.TransientSolution(2).HydrologySubTime,...
	      md.results.TransientSolution(3).SedimentHead,...
	      md.results.TransientSolution(3).EplHead,...
	      md.results.TransientSolution(3).SedimentHeadResidual};
