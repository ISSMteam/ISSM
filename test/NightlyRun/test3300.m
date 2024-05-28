md=triangle(model(),'../Exp/Square.exp',100000.);
md=setmask(md,'','');
md=parameterize(md,'../Par/SquareSheetConstrained.par');
md.transient=deactivateall(md.transient);
md.transient.ishydrology=1;
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',1);
md.hydrology=(hydrologydc);
md.hydrology=initialize(md.hydrology,md);
md.hydrology.isefficientlayer=1;
md.hydrology.sedimentlimit_flag=1;
md.hydrology.sedimentlimit=400.0;
md.hydrology.sediment_thickness=20.0;
md.initialization.sediment_head=0.0*ones(md.mesh.numberofvertices,1);
md.hydrology.spcsediment_head=NaN*ones(md.mesh.numberofvertices,1);
md.basalforcings.groundedice_melting_rate = 2.0*ones(md.mesh.numberofvertices,1);
md.hydrology.sediment_transmitivity=1.5e-4*ones(md.mesh.numberofvertices,1);

md.initialization.epl_head=0.0*ones(md.mesh.numberofvertices,1);
md.initialization.epl_thickness=1.0*ones(md.mesh.numberofvertices,1);
md.hydrology.spcepl_head=NaN*ones(md.mesh.numberofvertices,1);
md.hydrology.mask_eplactive_node=0*ones(md.mesh.numberofvertices,1);

md.hydrology.epl_conductivity=1.5e-2;
md.hydrology.epl_initial_thickness=1.0;
md.hydrology.epl_colapse_thickness=1.0e-6;
md.hydrology.epl_thick_comp=1;
md.hydrology.epl_max_thickness=5.0;

md.hydrology.transfer_flag=1.0;
md.hydrology.leakage_factor=3.9e-12;

times=0:0.002:8.0;
md.basalforcings.groundedice_melting_rate=ones(md.mesh.numberofvertices+1,length(times));

for i=1:length(times)
	if(times(i)<=1.0)
		md.basalforcings.groundedice_melting_rate(:,i)=1.0;
	elseif(times(i)<=6.0)
		md.basalforcings.groundedice_melting_rate(:,i)=-0.2;
	else
		md.basalforcings.groundedice_melting_rate(:,i)=0.0;
 end
end	

md.basalforcings.groundedice_melting_rate(end,:)=times;

md.timestepping.time_step=0.002;
md.timestepping.final_time=8.0;

md=solve(md,'Transient');

field_names     ={'SedimentWaterHead5','EplWaterHead5','SedimentWaterHead40','EplWaterHead40'};
field_tolerances={1e-13, 1e-13, 1e-13, 1e-13, 1e-13};
field_values={md.results.TransientSolution(5).SedimentHead, ...
							md.results.TransientSolution(5).EplHead,...
							md.results.TransientSolution(40).SedimentHead,...
							md.results.TransientSolution(40).EplHead};
