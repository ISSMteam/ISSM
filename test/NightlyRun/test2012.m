%Test Name: EarthSlc Hydrology loading

%mesh earth:
md=model;
load ../Data/SlcTestMesh.mat;
md.mesh=SlcMesh; %700 km resolution mesh

%Geometry for the bed, arbitrary thickness of 100: 
md.geometry.bed=zeros(md.mesh.numberofvertices,1);
md.geometry.base=md.geometry.bed;
md.geometry.thickness=100*ones(md.mesh.numberofvertices,1);
md.geometry.surface=md.geometry.bed+md.geometry.thickness;

%solidearth loading:  {{{
md.hydrology=hydrologytws();
md.hydrology.spcwatercolumn=[md.geometry.thickness;0];
md.smb.mass_balance=zeros(md.mesh.numberofvertices,1);

%antarctica
xe=md.mesh.x(md.mesh.elements)*[1;1;1]/3;
ye=md.mesh.y(md.mesh.elements)*[1;1;1]/3;
ze=md.mesh.z(md.mesh.elements)*[1;1;1]/3;
re=sqrt(xe.^2+ye.^2+ze.^2);

late=asind(ze./re);
longe=atan2d(ye,xe);
pos=find(late < -80);
md.hydrology.spcwatercolumn(md.mesh.elements(pos,:))=md.hydrology.spcwatercolumn(md.mesh.elements(pos,:))-100;
posant=pos;

%greenland
pos=find(late>60 & late<90 & longe>-75 & longe<-15);
md.hydrology.spcwatercolumn(md.mesh.elements(pos,:))=md.hydrology.spcwatercolumn(md.mesh.elements(pos,:))-100;
posgre=pos;

%elastic loading from love numbers:
md.solidearth.lovenumbers=lovenumbers('maxdeg',1000);

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

%masstransport:
md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
md.basalforcings.floatingice_melting_rate=zeros(md.mesh.numberofvertices,1);
md.initialization.vx=zeros(md.mesh.numberofvertices,1);
md.initialization.vy=zeros(md.mesh.numberofvertices,1);
md.initialization.sealevel=zeros(md.mesh.numberofvertices,1);
md.initialization.str=0;
md.initialization.watercolumn=zeros(md.mesh.numberofvertices,1);

%Materials: 
md.materials=materials('hydro');

%Miscellaneous
md.miscellaneous.name='test2012';

%Solution parameters
md.cluster.np=3;
md.solidearth.settings.reltol=NaN;
md.solidearth.settings.abstol=1e-3;
md.solidearth.settings.sealevelloading=1;
md.solidearth.settings.isgrd=1;
md.solidearth.settings.ocean_area_scaling=0;
md.solidearth.settings.grdmodel=1;

%Physics: 
md.transient.issmb=0; 
md.transient.isstressbalance=0;
md.transient.isthermal=0;
md.transient.ismasstransport=0;
md.transient.ishydrology=1;
md.transient.isslc=1;
md.solidearth.requested_outputs={'Sealevel','Bed'};

% max number of iteration reverted back to 10 (i.e., the original default value)
md.solidearth.settings.maxiter=10;

%eustatic run:
md.solidearth.settings.selfattraction=0;
md.solidearth.settings.elastic=0;
md.solidearth.settings.rotation=0;
md.solidearth.settings.viscous=0;

md=solve(md,'Transient');
Seustatic=md.results.TransientSolution.Sealevel;
Beustatic=md.results.TransientSolution.Bed;

%eustatic + selfattraction run:
md.solidearth.settings.selfattraction=1;
md.solidearth.settings.elastic=0;
md.solidearth.settings.rotation=0;
md.solidearth.settings.viscous=0;
md=solve(md,'tr');
Sselfattraction=md.results.TransientSolution.Sealevel;
Bselfattraction=md.results.TransientSolution.Bed;

%eustatic + selfattraction + elastic run:
md.solidearth.settings.selfattraction=1;
md.solidearth.settings.elastic=1;
md.solidearth.settings.rotation=0;
md.solidearth.settings.viscous=0;
md=solve(md,'tr');
Selastic=md.results.TransientSolution.Sealevel;
Belastic=md.results.TransientSolution.Bed;

%eustatic + selfattraction + elastic + rotation run:
md.solidearth.settings.selfattraction=1;
md.solidearth.settings.elastic=1;
md.solidearth.settings.rotation=1;
md.solidearth.settings.viscous=0;
md=solve(md,'tr');
Srotation=md.results.TransientSolution.Sealevel;
Brotation=md.results.TransientSolution.Bed;

%Fields and tolerances to track changes
field_names={'Seustatic','Sselfattraction','Selastic','Srotation','Beustatic','Bselfattraction','Belastic','Brotation'};
field_tolerances={1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13};
field_values={Seustatic,Sselfattraction,Selastic,Srotation,Beustatic,Bselfattraction,Belastic,Brotation};
