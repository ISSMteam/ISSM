%Test Name: RegionalMesh2dElasticGRDSpadaCap

%Flat projected mesh covering colatitudes 0--20 degrees.  The Spada load
%is a 10 degree polar cap, so the extra 10 degrees provide the response
%domain outside the load.
md=model;
re=md.solidearth.planetradius;
phi=(0:10:360)'*pi/180;

dom.x=cos(phi)*re*20*pi/180;
dom.y=sin(phi)*re*20*pi/180;
dom.x(end)=dom.x(1);
dom.y(end)=dom.y(1);
dom.nods=length(dom.x);

cap.x=cos(phi)*re*10*pi/180;
cap.y=sin(phi)*re*10*pi/180;
cap.x(end)=cap.x(1);
cap.y(end)=cap.y(1);
cap.nods=length(cap.x);
cap.closed=1;
cap.density=1;
cap.name='SpadaCap';
cap.Geometry='Polygon';
cap.BoundingBox=[min(cap.x) min(cap.y);max(cap.x) max(cap.y)];

md=bamg(md,'domain',dom,'subdomains',cap,'hmin',100e3,'hmax',100e3,...
	'KeepVertices',0,'MaxCornerAngle',1e-10,'NoBoundaryRefinement',1);

colat=sqrt(md.mesh.x.^2+md.mesh.y.^2)/re;
md.mesh.lat=90-colat*180/pi;
md.mesh.long=atan2d(md.mesh.y,md.mesh.x);

%Spada's 1500 m cap load, tapered to zero at 10 degrees colatitude.
load=zeros(md.mesh.numberofvertices,1);
pos=find(colat<=10*pi/180);
load(pos)=1500*sqrt((cos(colat(pos))-cosd(10))./(1-cosd(10)));

%Geometry and loading history.
md.geometry.bed=-zeros(md.mesh.numberofvertices,1);
md.geometry.base=md.geometry.bed;
md.geometry.thickness=0*ones(md.mesh.numberofvertices,1);
md.geometry.surface=md.geometry.bed+md.geometry.thickness;
md.masstransport.spcthickness=[md.geometry.thickness;0];
md.masstransport.spcthickness(1:end-1)=load;
md.smb.mass_balance=zeros(md.mesh.numberofvertices,1);

%Visco-elastic loading from the same temporal Love numbers as test2091.
load ../Data/lnb_temporal.mat
maxdeg=129;
mindeg=1;
md.solidearth.lovenumbers.h=ht;
md.solidearth.lovenumbers.h(maxdeg+1,:)=0;
md.solidearth.lovenumbers.h(1:mindeg,:)=0;
md.solidearth.lovenumbers.k=kt;
md.solidearth.lovenumbers.k(maxdeg+1,:)=-1;
md.solidearth.lovenumbers.k(1:mindeg,:)=-1;
md.solidearth.lovenumbers.l=lt;
md.solidearth.lovenumbers.l(maxdeg+1,:)=0;
md.solidearth.lovenumbers.l(1:mindeg,:)=0;
md.solidearth.lovenumbers.th=tht(1:maxdeg+1,:);
md.solidearth.lovenumbers.tk=tkt(1:maxdeg+1,:);
md.solidearth.lovenumbers.tl=tlt(1:maxdeg+1,:);
md.solidearth.lovenumbers.pmtf_colinear=pmtf1;
md.solidearth.lovenumbers.pmtf_ortho=pmtf2;
md.solidearth.lovenumbers.timefreq=time;

%Masks and regional deformation-only settings.
md.mask.ice_levelset=ones(md.mesh.numberofvertices,1);
md.mask.ice_levelset(pos)=-1;
md.mask.ocean_levelset=ones(md.mesh.numberofvertices,1);

md.timestepping.start_time=0;
md.timestepping.time_step=1000;
md.timestepping.final_time=12000;
time1=md.timestepping.start_time:md.timestepping.time_step:md.timestepping.final_time;
md.masstransport.spcthickness=repmat(md.masstransport.spcthickness,[1 length(time1)]);
md.masstransport.spcthickness(end,:)=time1;

md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
md.basalforcings.floatingice_melting_rate=zeros(md.mesh.numberofvertices,1);
md.initialization.vx=zeros(md.mesh.numberofvertices,1);
md.initialization.vy=zeros(md.mesh.numberofvertices,1);
md.initialization.sealevel=zeros(md.mesh.numberofvertices,1);
md.initialization.bottompressure=zeros(md.mesh.numberofvertices,1);
md.initialization.dsl=zeros(md.mesh.numberofvertices,1);
md.initialization.str=0;
md.materials=materials('hydro');
md.materials.rho_ice=931;

md.miscellaneous.name='test2095';
md.cluster=generic('name',oshostname(),'np',18);
md.solidearth.settings.isgrd=1;
md.solidearth.settings.grdmodel=1;
md.solidearth.settings.sealevelloading=0;
md.solidearth.settings.grdocean=0;
md.solidearth.settings.selfattraction=1;
md.solidearth.settings.elastic=1;
md.solidearth.settings.viscous=1;
md.solidearth.settings.rotation=0;
md.solidearth.settings.horiz=1;
md.solidearth.settings.ocean_area_scaling=0;
md.solidearth.settings.timeacc=md.timestepping.time_step;
md.solidearth.settings.degacc=.01;
md.solidearth.settings.viscoussampling=20;
md.solidearth.settings.maxiter=10;

md.transient.issmb=0;
md.transient.isstressbalance=0;
md.transient.isthermal=0;
md.transient.ismasstransport=1;
md.transient.isslc=1;
md.solidearth.requested_outputs={'SealevelGRD','BedGRD','BedNorthGRD','BedEastGRD',...
	'SealevelBarystaticIceLoad','SealevelBarystaticIceWeights','SealevelBarystaticIceArea',...
	'SealevelBarystaticIceMask','SealevelBarystaticIceLatbar','SealevelBarystaticIceLongbar'};
md.settings.results_on_nodes={'SealevelBarystaticIceWeights'};

md=solve(md,'Transient');

clear S B H E
for i=1:length(time1)-1
	S(:,i)=md.results.TransientSolution(i).SealevelGRD;
	B(:,i)=md.results.TransientSolution(i).BedGRD;
	H(:,i)=md.results.TransientSolution(i).BedNorthGRD;
	E(:,i)=md.results.TransientSolution(i).BedEastGRD;
end
S=cumsum(S,2);
B=cumsum(B,2);
H=cumsum(H,2);
E=cumsum(E,2);

field_names={'Bed','Sealevel','BedHorizontals'};
field_tolerances={1e-12,1e-12,1e-12};
field_values={B,S,H};
