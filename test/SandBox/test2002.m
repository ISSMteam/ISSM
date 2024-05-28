%Test Name: EarthSlc

%mesh earth:
md=model;
%md.mesh=gmshplanet('radius',6.371012*10^3,'resolution',700.); %700 km resolution mesh
load test2002_mesh 
md.mesh=mesh;

%Geometry for the bed, arbitrary thickness of 1000: 
md.geometry.bed=-ones(md.mesh.numberofvertices,1);
md.geometry.base=md.geometry.bed;
md.geometry.thickness=1000*ones(md.mesh.numberofvertices,1);
md.geometry.surface=md.geometry.bed+md.geometry.thickness;


%parameterize solidearth solution:
%solidearth loading:  {{{
md.masstransport.spcthickness=[md.geometry.thickness;0];
md.dsl.global_average_thermosteric_sea_level=[0;0];
md.dsl.sea_surface_height_above_geoid=zeros(md.mesh.numberofvertices+1,1);
md.dsl.sea_water_pressure_at_sea_floor=zeros(md.mesh.numberofvertices+1,1);
md.smb.mass_balance=zeros(md.mesh.numberofvertices,1);
%antarctica
late=sum(md.mesh.lat(md.mesh.elements),2)/3;
longe=sum(md.mesh.long(md.mesh.elements),2)/3;
%pos=find(late < -86);
pos=1136;
md.masstransport.spcthickness(md.mesh.elements(pos,:))= md.masstransport.spcthickness(md.mesh.elements(pos,:))-100;

vertices=md.mesh.elements(pos,:); 
pos1=find( md.mesh.elements(:,1)==vertices(1) | md.mesh.elements(:,2)==vertices(1) | md.mesh.elements(:,3)==vertices(1) );
pos2=find( md.mesh.elements(:,1)==vertices(2) | md.mesh.elements(:,2)==vertices(2) | md.mesh.elements(:,3)==vertices(2) );
pos3=find( md.mesh.elements(:,1)==vertices(3) | md.mesh.elements(:,2)==vertices(3) | md.mesh.elements(:,3)==vertices(3) );
polarelements=unique([pos1;pos2;pos3]);


%greenland
%pos=find(late>70 & late<80 & longe>-60 & longe<-30);
%md.masstransport.spcthickness(md.mesh.elements(pos,:))= md.masstransport.spcthickness(md.mesh.elements(pos,:))-100;

%elastic loading from love numbers:
md.solidearth.lovenumbers=lovenumbers('maxdeg',100);

%}}}
%mask:  {{{
md.mask.ice_levelset=ones(md.mesh.numberofvertices,1);
md.mask.ocean_levelset=-ones(md.mesh.numberofvertices,1);
pos=find(md.mesh.lat<0);
md.mask.ice_levelset(pos)=-1;
md.mask.ocean_levelset(pos)=1;

md.mask.ocean_levelset(md.mesh.elements(polarelements,:))=-1;
md.mask.ocean_levelset(md.mesh.elements(1136,:))=1;

md.mask.ice_levelset(md.mesh.elements(polarelements,:))=1;
md.mask.ice_levelset(md.mesh.elements(1136,:))=-1;


% }}}

%time stepping: 
md.timestepping.start_time=0;
md.timestepping.time_step=1;
md.timestepping.final_time=1;


md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
md.basalforcings.floatingice_melting_rate=zeros(md.mesh.numberofvertices,1);
md.initialization.vx=zeros(md.mesh.numberofvertices,1);
md.initialization.vy=zeros(md.mesh.numberofvertices,1);
md.initialization.sealevel=zeros(md.mesh.numberofvertices,1);
md.initialization.bottompressure=zeros(md.mesh.numberofvertices,1);
md.initialization.dsl=zeros(md.mesh.numberofvertices,1);
md.initialization.str=0;

%Materials: 
md.materials=materials('hydro');


%Miscellaneous
md.miscellaneous.name='test2002';

%Solution parameters
md.solidearth.settings.reltol=NaN;
md.solidearth.settings.abstol=1e-3;
md.solidearth.settings.computesealevelchange=1;
md.solidearth.settings.isgrd=1;
md.solidearth.settings.ocean_area_scaling=0;
md.solidearth.settings.grdmodel=1;

%Physics: 
md.transient.issmb=0; 
md.transient.isstressbalance=0;
md.transient.isthermal=0;
md.transient.ismasstransport=1;
md.transient.isoceantransport=1;
md.transient.isslc=1;
md.solidearth.requested_outputs={'Sealevel','SealevelBarystaticMask','DeltaDsl','DeltaBottomPressure','DeltaStr','DeltaIceThickness'};


% max number of iteration reverted back to 10 (i.e., the original default value)
md.solidearth.settings.maxiter=10;

%eustatic run:
md.solidearth.settings.rigid=0;
md.solidearth.settings.elastic=0;
md.solidearth.settings.rotation=0;
md=solve(md,'Transient');
Seustatic=md.results.TransientSolution.Sealevel;

%eustatic + rigid run:
md.solidearth.settings.rigid=1;
md.solidearth.settings.elastic=0;
md.solidearth.settings.rotation=0;
md=solve(md,'tr');
Srigid=md.results.TransientSolution.Sealevel;

%eustatic + rigid + elastic run:
md.solidearth.settings.rigid=1;
md.solidearth.settings.elastic=1;
md.solidearth.settings.rotation=0;

md.solidearth.requested_outputs={'Sealevel','SealevelBarystaticMask','DeltaDsl','DeltaBottomPressure','DeltaStr','DeltaIceThickness','SealevelchangeG','SealevelBarystaticOceanMask'};
md=solve(md,'tr');
Selastic=md.results.TransientSolution.Sealevel;

Seustatic_new=Seustatic;
Srigid_new=Srigid;
Selastic_new=Selastic;
save ../../../trunk-clean/test/NightlyRun/test2002_results Seustatic_new Srigid_new Selastic_new

%eustatic + rigid + elastic + rotation run:
md.solidearth.settings.rigid=1;
md.solidearth.settings.elastic=1;
md.solidearth.settings.rotation=1;
md=solve(md,'tr');
Srotation=md.results.TransientSolution.Sealevel;
error;

%Fields and tolerances to track changes
field_names={'Eustatic','Rigid','Elastic','Rotation'};
field_tolerances={1e-13,1e-13,1e-13,1e-13};
field_values={Seustatic,Srigid,Selastic,Srotation};
