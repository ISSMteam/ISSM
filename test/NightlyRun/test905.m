%Test Name: SquareNoDynHydrologyDCSmbCoupled
md=triangle(model(),'../Exp/Square.exp',100000.);
md=setmask(md,'','');
md=parameterize(md,'../Par/SquareNoDyn.par');
md.cluster = generic('name',oshostname(),'np',1);

md.transient.ishydrology=1;
md.transient.issmb=1;
md.hydrology=(hydrologydc);
md.hydrology=initialize(md.hydrology,md);
md.smb=(SMBgradientscomponents);

md.hydrology.isefficientlayer=1;

md.hydrology.sedimentlimit_flag=1;
md.hydrology.sedimentlimit=400.0;
md.hydrology.sediment_transmitivity=3.0*ones(md.mesh.numberofvertices,1);
md.hydrology.mask_thawed_node=ones(md.mesh.numberofvertices,1);

md.hydrology.mask_eplactive_node=zeros(md.mesh.numberofvertices,1);
md.hydrology.epl_conductivity=3.;
md.hydrology.epl_initial_thickness=20;
md.hydrology.epl_colapse_thickness=1.0e-3;
md.hydrology.epl_thick_comp=0;
md.hydrology.epl_max_thickness=1;

md.hydrology.spcsediment_head=NaN*ones(md.mesh.numberofvertices,1);
md.hydrology.spcepl_head=NaN*ones(md.mesh.numberofvertices,1);

md.initialization.sediment_head=zeros(md.mesh.numberofvertices,1);
md.initialization.epl_head=zeros(md.mesh.numberofvertices,1);
md.initialization.epl_thickness=ones(md.mesh.numberofvertices,1);

md.hydrology.steps_per_step=5;
md.smb.steps_per_step=10;
md.timestepping.time_step=1.;
md.timestepping.final_time=20.0;

smb_step=md.timestepping.time_step/md.smb.steps_per_step;
duration=[md.timestepping.start_time:smb_step:md.timestepping.final_time];

ddf=10.0e-3;
md.smb.accuref=[[0.5 0.5];[md.timestepping.start_time md.timestepping.final_time]];
md.smb.accualti=0.0;
md.smb.accugrad=[[0. 0.];[md.timestepping.start_time md.timestepping.final_time]];

md.smb.runoffref=0.9*duration*ddf;
md.smb.runoffref=[md.smb.runoffref;duration];
md.smb.runoffalti=0.0;
md.smb.runoffgrad=[[-6.5e-3*ddf -6.5e-3*ddf];[md.timestepping.start_time md.timestepping.final_time]];  %lapse rate *ddf*day per year

md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
md.basalforcings.floatingice_melting_rate=zeros(md.mesh.numberofvertices,1);

md=solve(md,'Transient');

field_names={'SedimentWaterHead1','EplWaterHead1','SedimentHeadResidual1',...
	     'SedimentWaterHead4','EplWaterHead4','SedimentHeadResidual4',...
	     'SedimentWaterHead5','EplWaterHead5','SedimentHeadResidual5',...
	     'SedimentWaterHead9','EplWaterHead9','SedimentHeadResidual9',...
	     'EplWaterHead20','EplWaterHeadSubstep20','SedimentWaterHead20',...
	     'SedimentWaterHeadSubstep20'};
field_tolerances={1e-13,1e-13,1e-13,...
		  1e-13,1e-13,1e-13,...
		  1e-13,5e-12,1e-11,...
		  1e-13,5e-12,1e-11,...
		  1e-13,1e-13,1e-13,...
		  1e-13};
field_values={md.results.TransientSolution(1).SedimentHead,...
	      md.results.TransientSolution(1).EplHead,...
	      md.results.TransientSolution(1).SedimentHeadResidual,...
	      md.results.TransientSolution(4).SedimentHead,...
	      md.results.TransientSolution(4).EplHead,...
	      md.results.TransientSolution(4).SedimentHeadResidual,...
	      md.results.TransientSolution(5).SedimentHead,...
	      md.results.TransientSolution(5).EplHead,...
	      md.results.TransientSolution(5).SedimentHeadResidual,...
	      md.results.TransientSolution(9).SedimentHead,...
	      md.results.TransientSolution(9).EplHead,...
	      md.results.TransientSolution(9).SedimentHeadResidual,...
	      md.results.TransientSolution(end).EplHead,...
	      md.results.TransientSolution(end).EplHeadSubstep,...
	      md.results.TransientSolution(end).SedimentHead,...
	      md.results.TransientSolution(end).SedimentHeadSubstep};
