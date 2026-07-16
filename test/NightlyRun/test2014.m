%Test Name: SlcBarystaticGroundinglineMigration

%mesh earth:
md=model;
load ../Data/SlcTestMesh.mat;
md.mesh=SlcMesh; %700 km resolution mesh

%Create a marine ice patch that starts grounded and becomes floating.
H0=1300;
H1=900;
bed=-1000;

md.geometry.bed=bed*ones(md.mesh.numberofvertices,1);
md.geometry.base=md.geometry.bed;
md.geometry.thickness=H0*ones(md.mesh.numberofvertices,1);
md.geometry.surface=md.geometry.base+md.geometry.thickness;

md.smb.mass_balance=zeros(md.mesh.numberofvertices,1);
md.dsl.global_average_thermosteric_sea_level=[0;0];
md.dsl.sea_surface_height_above_geoid=zeros(md.mesh.numberofvertices+1,1);
md.dsl.sea_water_pressure_at_sea_floor=zeros(md.mesh.numberofvertices+1,1);

xe=md.mesh.x(md.mesh.elements)*[1;1;1]/3;
ye=md.mesh.y(md.mesh.elements)*[1;1;1]/3;
ze=md.mesh.z(md.mesh.elements)*[1;1;1]/3;
re=sqrt(xe.^2+ye.^2+ze.^2);

late=asind(ze./re);
longe=atan2d(ye,xe);
posice=find(late>60 & late<90 & longe>-75 & longe<-15);

%mask:
mask=gmtmask(md.mesh.lat,md.mesh.long);
oceanmask=-ones(md.mesh.numberofvertices,1);
pos=find(mask==0); oceanmask(pos)=1;

icemask=ones(md.mesh.numberofvertices,1);
icemask(md.mesh.elements(posice,:))=-1;
oceanmask(md.mesh.elements(posice,:))=1;

md.mask.ice_levelset=icemask;
md.mask.ocean_levelset=oceanmask;

%prescribed thickness: first column is the initial state, second column
%drops the marine patch below flotation.
Hfinal=md.geometry.thickness;
Hfinal(icemask<0)=H1;
md.masstransport.spcthickness=[md.geometry.thickness Hfinal;0 1];

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
md.miscellaneous.name='test2014';

%Solution parameters
md.cluster.np=3;
md.solidearth.lovenumbers=lovenumbers('maxdeg',100);
md.solidearth.settings.reltol=NaN;
md.solidearth.settings.abstol=1e-3;
md.solidearth.settings.sealevelloading=1;
md.solidearth.settings.isgrd=1;
md.solidearth.settings.ocean_area_scaling=0;
md.solidearth.settings.grdmodel=1;
md.solidearth.settings.maxiter=10;

%Physics:
md.transient.issmb=0;
md.transient.isstressbalance=0;
md.transient.isthermal=0;
md.transient.ismasstransport=1;
md.transient.isgroundingline=1;
md.transient.isslc=1;
md.groundingline.migration='SubelementMigration';
md.solidearth.requested_outputs={'Sealevel','Bed','Thickness','MaskOceanLevelset','MaskIceLevelset',...
	'SealevelBarystaticIceLoad','SealevelBarystaticSeaLevelLoad',...
	'SealevelBarystaticIceWeights','SealevelBarystaticIceArea',...
	'SealevelBarystaticOceanWeights','SealevelBarystaticOceanArea','SealevelBarystaticOceanMask'};
md.settings.results_on_nodes={'SealevelBarystaticIceWeights','SealevelBarystaticOceanWeights'};

%Keep this test focused on barystatic load accounting.
md.solidearth.settings.selfattraction=0;
md.solidearth.settings.elastic=0;
md.solidearth.settings.rotation=0;
md.solidearth.settings.viscous=0;
md=solve(md,'tr');

results=md.results.TransientSolution(1);
rhoi=md.materials.rho_ice;
rhow=md.materials.rho_water;

oceanarea=sum(results.SealevelBarystaticOceanArea);
icearea=results.SealevelBarystaticIceArea;
iceweights=results.SealevelBarystaticIceWeights;

H0v=md.geometry.thickness(md.mesh.elements);
bed0=md.geometry.bed(md.mesh.elements);
sl0=md.initialization.sealevel(md.mesh.elements);
waterdepth0=max(sl0-bed0,0);
haf0=H0v-rhow/rhoi*waterdepth0;

H0avg=sum(H0v.*iceweights,2);
haf0avg=sum(haf0.*iceweights,2);

iceload=sum(results.SealevelBarystaticIceLoad);
sealevelload=sum(results.SealevelBarystaticSeaLevelLoad);
netload=iceload+sealevelload;

iceload_geom=-sum(H0avg.*icearea)*rhoi;
sealevelload_geom=sum((H0avg-haf0avg).*icearea)*rhoi;
netload_geom=-sum(haf0avg.*icearea)*rhoi;

bslcice=results.BslcIce;
bslcice_geom=-iceload_geom/rhow/oceanarea;
bslcice_load=-iceload/rhow/oceanarea;

bslcsealevelload=-sealevelload/rhow/oceanarea;
bslcsealevelload_geom=-sealevelload_geom/rhow/oceanarea;

bslcnetload=-netload/rhow/oceanarea;
bslcnetload_geom=-netload_geom/rhow/oceanarea;

bslcice_diff=single(bslcice_geom)-single(bslcice);
bslcice_load_diff=single(bslcice_load)-single(bslcice);
bslcsealevelload_diff=single(bslcsealevelload_geom)-single(bslcsealevelload);
bslcnetload_diff=single(bslcnetload_geom)-single(bslcnetload);
cumbslcice_diff=single(results.CumBslcIce)-single(bslcice);

%Fields and tolerances to track changes
field_names={'BslcIce','BslcIceGeom','BslcIceDiff','BslcIceLoadDiff',...
	'BslcSeaLevelLoad','BslcSeaLevelLoadGeom','BslcSeaLevelLoadDiff',...
	'BslcNetLoad','BslcNetLoadGeom','BslcNetLoadDiff','CumBslcIceDiff'};
field_tolerances={2e-12,2e-12,1e-13,1e-13,...
	2e-12,2e-12,1e-13,...
	2e-12,2e-12,1e-13,1e-13};
field_values={bslcice,bslcice_geom,bslcice_diff,bslcice_load_diff,...
	bslcsealevelload,bslcsealevelload_geom,bslcsealevelload_diff,...
	bslcnetload,bslcnetload_geom,bslcnetload_diff,cumbslcice_diff};
