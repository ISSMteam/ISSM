%Test Name: RoundSheetShelfGLMigrationMOLHO2d
radius=1.e6;
shelfextent=2.e5;

md=roundmesh(model(),radius,50000.);
%fix center node to 0,0
rad=sqrt(md.mesh.x.^2+md.mesh.y.^2);
pos=find(rad==min(rad));
md.mesh.x(pos)=0.; md.mesh.y(pos)=0.; %the closest node to the center is changed to be exactly at the center
xelem=mean(md.mesh.x(md.mesh.elements),2);
yelem=mean(md.mesh.y(md.mesh.elements),2);
rad=sqrt(xelem.^2+yelem.^2);
flags=zeros(md.mesh.numberofelements,1);
pos=find(rad>=(radius-shelfextent));
flags(pos)=1;
md=setmask(md,flags,''); 
md=parameterize(md,'../Par/RoundSheetShelf.par');
md=setflowequation(md,'MOLHO','all');
md.cluster=generic('name',oshostname(),'np',3);

md.transient.isthermal=0;
md.transient.ismasstransport=0;
md.transient.issmb=1;
md.transient.isstressbalance=1;
md.transient.isgroundingline=1;

md=SetMOLHOBC(md);
%test different grounding line dynamics.
md.groundingline.migration='AggressiveMigration';
md=solve(md,'Transient');
element_on_iceshelf_agressive=(md.results.TransientSolution.MaskOceanLevelset);
vel_agressive=(md.results.TransientSolution.Vel);

md.groundingline.migration='SoftMigration';
md=solve(md,'Transient');
element_on_iceshelf_soft=(md.results.TransientSolution.MaskOceanLevelset);
vel_soft=(md.results.TransientSolution.Vel);

md.mask.ocean_levelset=md.geometry.thickness + md.materials.rho_water/md.materials.rho_ice*md.geometry.bed;
md.groundingline.migration='SubelementMigration';
md.groundingline.friction_interpolation='SubelementFriction1';
md=solve(md,'Transient');
element_on_iceshelf_subelement=(md.results.TransientSolution.MaskOceanLevelset);
vel_subelement=(md.results.TransientSolution.Vel);

md.groundingline.friction_interpolation='SubelementFriction2';
md=solve(md,'Transient');
element_on_iceshelf_subelement2=(md.results.TransientSolution.MaskOceanLevelset);
vel_subelement2=(md.results.TransientSolution.Vel);

%Fields and tolerances to track changes
field_names     ={'ElementOnIceShelfAggressive','VelAggressive','ElementOnIceShelfSoft','VelSoft','ElementOnIceShelfSubelement','VelSubelement','ElementOnIceShelfSubelement2','VelSubelement2'};
field_tolerances={1e-13,2e-12,1e-13,2e-12,1e-13,3e-12,1e-13,2e-12};
field_values={element_on_iceshelf_agressive,vel_agressive,element_on_iceshelf_soft,vel_soft,element_on_iceshelf_subelement,vel_subelement,element_on_iceshelf_subelement2,vel_subelement2};
