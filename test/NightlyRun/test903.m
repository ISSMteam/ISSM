%Test Name: SquareNoDynExtrudedHydrologyDCOneLayer
md=triangle(model(),'../Exp/Square.exp',100000.);
md=setmask(md,'','');
md=parameterize(md,'../Par/SquareNoDyn.par');
md.cluster=generic('name',oshostname(),'np',1);

md.transient.ishydrology=1;
md.hydrology=(hydrologydc);
md.hydrology=initialize(md.hydrology,md);

md.hydrology.isefficientlayer=0;
md.hydrology.sedimentlimit_flag=1;
md.hydrology.sedimentlimit=8000.0;
md.hydrology.mask_thawed_node=ones(md.mesh.numberofvertices,1);
md.initialization.sediment_head=0.0*ones(md.mesh.numberofvertices,1);
md.hydrology.spcsediment_head=NaN*ones(md.mesh.numberofvertices,1);
pos=find(md.mesh.y==0);
md.hydrology.spcsediment_head(pos)=0.0;

md.basalforcings.groundedice_melting_rate = 2.0*ones(md.mesh.numberofvertices,1);
md.basalforcings.floatingice_melting_rate = 0.0*ones(md.mesh.numberofvertices,1);
md.hydrology.sediment_transmitivity= 3.0*ones(md.mesh.numberofvertices,1);

md.timestepping.time_step=0;
md.timestepping.final_time=1.0;
md=extrude(md,3,1.1);
md=solve(md,'Hydrology');

%Fields and tolerances to track changes
%you can also compare with an analitic solution, but it is exact
%only if no limits are applied
%analitic=(md.mesh.y.^2-2*md.mesh.y*1.0e6)*(-2.0/(2*md.constants.yts*md.hydrology.sediment_transmitivity))
field_names     ={'SedimentWaterHead','SedimentHeadResidual'};
field_tolerances={1e-13, 3e-10};
field_values={md.results.HydrologySolution.SedimentHead,md.results.HydrologySolution.SedimentHeadResidual};
