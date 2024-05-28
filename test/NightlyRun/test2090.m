%Test Name: ViscoElasticEarthSlc

%mesh earth:
md=model;
load ../Data/SlcTestMesh.mat;
md.mesh=SlcMesh; %700 km resolution mesh

%Geometry for the bed, arbitrary thickness of 1000: 
md.geometry.bed=-ones(md.mesh.numberofvertices,1);
md.geometry.base=md.geometry.bed;
md.geometry.thickness=100*ones(md.mesh.numberofvertices,1);
md.geometry.surface=md.geometry.bed+md.geometry.thickness;


%parameterize solidearth solution:
%solidearth loading:  {{{
md.masstransport.spcthickness=[md.geometry.thickness;0];

md.smb.mass_balance=zeros(md.mesh.numberofvertices,1);
%antarctica
xe=md.mesh.x(md.mesh.elements)*[1;1;1]/3;
ye=md.mesh.y(md.mesh.elements)*[1;1;1]/3;
ze=md.mesh.z(md.mesh.elements)*[1;1;1]/3;
re=sqrt(xe.^2+ye.^2+ze.^2);

late=asind(ze./re);
longe=atan2d(ye,xe);
pos=find(late < -80);
%pos=10;
md.masstransport.spcthickness(md.mesh.elements(pos,:))=md.masstransport.spcthickness(md.mesh.elements(pos,:)) -100;
posant=pos;
%greenland
pos=find(late>60 & late<90 & longe>-75 & longe<-15);
md.masstransport.spcthickness(md.mesh.elements(pos,:))=md.masstransport.spcthickness(md.mesh.elements(pos,:)) -100;
posgre=pos;

%visco-elastic loading from love numbers that we load ourselves. We 
%still use the lovenumbers constructor to initialize the fields: 
load ../Data/lnb_temporal.mat
%md.solidearth.lovenumbers=lovenumbers('maxdeg',1000);
%md.solidearth.lovenumbers.timefreq=[0];

md.solidearth.lovenumbers.h=[ht];
md.solidearth.lovenumbers.k=[kt];
md.solidearth.lovenumbers.l=[lt];
%md.solidearth.lovenumbers.h=repmat(md.solidearth.lovenumbers.h,1,100);
%md.solidearth.lovenumbers.k=repmat(md.solidearth.lovenumbers.k,1,100);
%md.solidearth.lovenumbers.l=repmat(md.solidearth.lovenumbers.l,1,100);
md.solidearth.lovenumbers.th=tht;
md.solidearth.lovenumbers.tk=tkt;
md.solidearth.lovenumbers.tl=tlt;
md.solidearth.lovenumbers.pmtf_colinear=pmtf1;
md.solidearth.lovenumbers.pmtf_ortho=pmtf2;
md.solidearth.lovenumbers.timefreq=time;

%}}}
%mask:  {{{
mask=gmtmask(md.mesh.lat,md.mesh.long);
oceanmask=-ones(md.mesh.numberofvertices,1);
pos=find(mask==0); oceanmask(pos)=1;

icemask=ones(md.mesh.numberofvertices,1);
icemask(md.mesh.elements(posant))=-1;
icemask(md.mesh.elements(posgre))=-1;
%oceanmask(md.mesh.elements(posant))=1;

md.mask.ice_levelset=icemask;
md.mask.ocean_levelset=oceanmask;
% }}}

%time stepping: 
md.timestepping.start_time=0;
md.timestepping.time_step=100;
md.timestepping.final_time=1000;

time1=md.timestepping.start_time:md.timestepping.time_step:md.timestepping.final_time;
md.masstransport.spcthickness=repmat(md.masstransport.spcthickness, [1 length(time1)]);
md.masstransport.spcthickness(1:end-1,3:end)=0;
md.masstransport.spcthickness(end,:)=time1;

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
md.miscellaneous.name='test2090';

%Solution parameters
md.cluster.np=10;
md.solidearth.settings.reltol=NaN;
md.solidearth.settings.abstol=1e-3;
md.solidearth.settings.sealevelloading=1;
md.solidearth.settings.isgrd=1;
md.solidearth.settings.ocean_area_scaling=0;
md.solidearth.settings.grdmodel=1;
md.solidearth.settings.timeacc=md.timestepping.time_step;
md.solidearth.settings.degacc=.1;

%Physics:
md.transient.issmb=0; 
md.transient.isstressbalance=0;
md.transient.isthermal=0;
md.transient.ismasstransport=1;
md.transient.isslc=1;
md.solidearth.requested_outputs={'Sealevel','Bed'};

% max number of iteration reverted back to 10 (i.e., the original default value)
md.solidearth.settings.maxiter=10;

%eustatic + selfattraction + elastic + rotation run:
md.solidearth.settings.selfattraction=1;
md.solidearth.settings.elastic=1;
md.solidearth.settings.viscous=1;
md.solidearth.settings.rotation=1;
md=solve(md,'tr');


clear S B
for i=1:length(time1)-1
S(:,i)=md.results.TransientSolution(i).Sealevel;
B(:,i)=md.results.TransientSolution(i).Bed;
end

hold on
plot(md.timestepping.start_time:md.timestepping.time_step:(md.timestepping.final_time-md.timestepping.time_step), [B(916,:)]);

%Fields and tolerances to track changes
field_names={'Sealevel','Bed'};
field_tolerances={1e-13,1e-13};
field_values={S,B};
