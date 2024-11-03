%Test Name: SlcBarystatic

%mesh earth:
md=model;
load ../Data/SlcTestMesh.mat;
md.mesh=SlcMesh; %700 km resolution mesh

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
xe=md.mesh.x(md.mesh.elements)*[1;1;1]/3;
ye=md.mesh.y(md.mesh.elements)*[1;1;1]/3;
ze=md.mesh.z(md.mesh.elements)*[1;1;1]/3;
re=sqrt(xe.^2+ye.^2+ze.^2);

late=asind(ze./re);
longe=atan2d(ye,xe);
pos=find(late < -80);
md.masstransport.spcthickness(md.mesh.elements(pos,:))= md.masstransport.spcthickness(md.mesh.elements(pos,:))-100;
posant=pos;
%greenland
pos=find(late>60 & late<90 & longe>-75 & longe<-15);
md.masstransport.spcthickness(md.mesh.elements(pos,:))= md.masstransport.spcthickness(md.mesh.elements(pos,:))-100;
posgre=pos;

%elastic loading from love numbers:
md.solidearth.lovenumbers=lovenumbers('maxdeg',100);

%}}}
%mask:  {{{
mask=gmtmask(md.mesh.lat,md.mesh.long);
oceanmask=-ones(md.mesh.numberofvertices,1);
pos=find(mask==0); oceanmask(pos)=1;

icemask=ones(md.mesh.numberofvertices,1);
icemask(md.mesh.elements(posant))=-1;
icemask(md.mesh.elements(posgre))=-1;

md.mask.ice_levelset=icemask;
md.mask.ocean_levelset=oceanmask;
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
md.miscellaneous.name='test2011';

%Solution parameters
md.cluster.np=3;
md.solidearth.settings.reltol=1e-10;
md.solidearth.settings.abstol=NaN; 
md.solidearth.settings.sealevelloading=1;
md.solidearth.settings.isgrd=1;
md.solidearth.settings.ocean_area_scaling=0;
md.solidearth.settings.grdmodel=1;

%Physics: 
md.transient.issmb=0; 
md.transient.isstressbalance=0;
md.transient.isthermal=0;
md.transient.ismasstransport=1;
md.transient.isslc=1;
md.solidearth.requested_outputs={'Sealevel','DeltaIceThickness','SealevelBarystaticIceMask','SealevelBarystaticHydroMask','SealevelBarystaticBpMask','Bed','SealevelBarystaticIceWeights','SealevelBarystaticOceanWeights','SealevelBarystaticIceArea','SealevelBarystaticOceanArea','SealevelBarystaticOceanMask','SealevelGRD','SealevelBarystaticOceanLongbar'};
md.settings.results_on_nodes={'SealevelBarystaticIceWeights','SealevelBarystaticOceanWeights'};

%eustatic + rigid + elastic + rotation run:
md.solidearth.settings.selfattraction=1;
md.solidearth.settings.elastic=1;
md.solidearth.settings.rotation=0;
md.solidearth.settings.viscous=0;
md=solve(md,'tr');

%recover barystatic: 
results=md.results.TransientSolution(1); 
bslc1=results.BslcIce; 

%alternative way of computing barystatic: 
icearea=results.SealevelBarystaticIceArea;
iceweights=results.SealevelBarystaticIceWeights; 
dH=results.DeltaIceThickness(md.mesh.elements); 
dHavg=sum(dH.*iceweights,2);
totaloceanarea=sum(results.SealevelBarystaticOceanArea);
bslc2=-sum(dHavg.*icearea)*md.materials.rho_ice/md.materials.rho_water/totaloceanarea;

%another way of computing barystatic:
oceanarea=results.SealevelBarystaticOceanArea; 
oceanweights=results.SealevelBarystaticOceanWeights; 
rsl=results.Sealevel-(results.Bed-md.geometry.bed); 
rsl=rsl(md.mesh.elements); 
rslavg=sum(rsl.*oceanweights,2);
bslc3=sum(rslavg.*oceanarea)/totaloceanarea; 

%need to change precision before subtraction because of differences in results 
%at high precision under macOS versus Linux (print values of bslc and bslc2 to %verify that the difference is negligible)
bslc_diff21=single(bslc2)-single(bslc1);
bslc_diff31=single(bslc3)-single(bslc2);

%Fields and tolerances to track changes
field_names={'BarystaticIce1','BarystaticIce2','BarystaticIceDiff21','BarystaticIce3','BarystaticIceDiff31'};
field_tolerances={1e-13,1e-13,1e-13,1e-13,1e-13};
field_values={bslc1,bslc2,bslc_diff21,bslc3,bslc_diff31};
