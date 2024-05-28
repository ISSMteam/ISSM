%Test Name: PolarMotion

%mesh earth:
md=model;md.mesh=mesh3dsurface();
re= md.solidearth.planetradius;
phi1=(0:10:360)/180*pi;
phi2=(0:10:360)/180*pi;

dom=struct;
dom.x=cos(phi1')*re*pi;
dom.y=sin(phi1')*re*pi;dom.y(end)=dom.y(1);
dom.nods=length(dom.x);   

co=struct;
rad=re*pi*[linspace(1e-1,1,3)/18];
for i=1:length(rad);
co(i).x=cos(phi2')*rad(i);
co(i).y=sin(phi2')*rad(i);co(i).y(end)=co(i).y(1);
co(i).nods=length(co(i).x);co(i).closed=1;co(i).density=1;co(i).name='Icedisklimit';co(i).Geometry='Polygon';
co(i).BoundingBox=[min(co(i).x) min(co(i).y); max(co(i).x) max(co(i).y)];
end

defaultoptions={'KeepVertices',0,'MaxCornerAngle',0.0000000001,'NoBoundaryRefinement',1}; 
md2d=bamg(model,'domain',dom,'subdomains',co,'hmin',50e3,'hmax',10000e3,defaultoptions{:});

colat=sqrt(md2d.mesh.x.^2+md2d.mesh.y.^2)/re;
lon=atan2(md2d.mesh.y,md2d.mesh.x);

x=sin(colat).*cos(lon)*re;
y=sin(colat).*sin(lon)*re;
z=cos(colat)*re;

%load
longL=75;
latL=65;

%rotate mesh around load center
R1=[cosd(90-latL) 0 sind(90-latL) ;0 1 0; -sind(90-latL) 0 cosd(90-latL)];
R2=[cosd(longL) -sind(longL) 0; sind(longL) cosd(longL) 0; 0 0 1];

for i=1:length(x);
	coord= R2*R1*[x(i);y(i);z(i)];
	x(i)=coord(1);
	y(i)=coord(2);
	z(i)=coord(3);
end

colat=acos(z/re);
lon=atan2(y,x);


md.mesh.x=x;md.mesh.y=y;md.mesh.z=z;md.mesh.lat=90-colat*180/pi;md.mesh.long=lon*180/pi;
md.mesh.r=sqrt(x.^2+y.^2+z.^2);
md.mesh.numberofelements=md2d.mesh.numberofelements;md.mesh.numberofvertices=md2d.mesh.numberofvertices;
md.mesh.elements=md2d.mesh.elements;
md.mesh.vertexconnectivity=NodeConnectivity(md.mesh.elements,md.mesh.numberofvertices);
%md.mesh.area=averaging(md,GetAreas3DTria(md.mesh.elements,md.mesh.x,md.mesh.y,md.mesh.z),4);
md.mesh.area=GetAreas3DTria(md.mesh.elements,md.mesh.x,md.mesh.y,md.mesh.z);

%Geometry for the bed
md.geometry.bed=-zeros(md.mesh.numberofvertices,1);
md.geometry.base=md.geometry.bed;
md.geometry.thickness=0*ones(md.mesh.numberofvertices,1);
md.geometry.surface=md.geometry.bed+md.geometry.thickness;


%parameterize solidearth solution:
%solidearth loading:  {{{
md.masstransport.spcthickness=[md.geometry.thickness;0];
md.materials.rho_ice=931;


md.smb.mass_balance=zeros(md.mesh.numberofvertices,1);
%antarctica
xe=md.mesh.x(md.mesh.elements)*[1;1;1]/3;
ye=md.mesh.y(md.mesh.elements)*[1;1;1]/3;
ze=md.mesh.z(md.mesh.elements)*[1;1;1]/3;
re=sqrt(xe.^2+ye.^2+ze.^2);

late=asind(ze./re);
longe=atan2d(ye,xe);


longi=md.mesh.long;
lati=md.mesh.lat;
delPhi=abs(lati-latL); delLambda=abs(longi-longL); if (delLambda>pi)delLambda=2*pi-delLambda; end
alpha=2.*asin(sqrt(sind(delPhi/2).^2+cosd(lati).*cosd(latL).*sind(delLambda/2).^2));

pos=find(alpha<=10/180*pi-0.01);

%mass=3.607171340900778E+018;
area_element=GetAreasSphericalTria(md.mesh.elements,md.mesh.lat,md.mesh.long,md.solidearth.planetradius);

%md.masstransport.spcthickness(md.mesh.elements(pos,:))=mass/md.materials.rho_ice/area_element(pos);
md.masstransport.spcthickness(pos,1)=1500 * sqrt( (cos(alpha(pos)) - cosd(10)) /(1-cosd(10)) );


%visco-elastic loading from love numbers that we load ourselves. We 
%still use the lovenumbers constructor to initialize the fields: 
load ../Data/lnb_temporal.mat
%md.solidearth.lovenumbers=lovenumbers('maxdeg',1000);
%md.solidearth.lovenumbers.timefreq=[0];

maxdeg=3;
mindeg=2;

md.solidearth.lovenumbers.h=ht;
md.solidearth.lovenumbers.h(maxdeg+1,:)=0;
%md.solidearth.lovenumbers.h(maxdeg+2,:)=0;
%md.solidearth.lovenumbers.h(maxdeg+3:end,:)=[];
md.solidearth.lovenumbers.h(1:mindeg,:)=0;
md.solidearth.lovenumbers.k=kt;
md.solidearth.lovenumbers.k(maxdeg+1,:)=-1;
%md.solidearth.lovenumbers.k(maxdeg+2,:)=0;
%md.solidearth.lovenumbers.k(maxdeg+3:end,:)=[];
md.solidearth.lovenumbers.k(1:mindeg,:)=-1;
md.solidearth.lovenumbers.l=lt;
md.solidearth.lovenumbers.l(maxdeg+1,:)=0;
%md.solidearth.lovenumbers.l(maxdeg+2,:)=0;
%md.solidearth.lovenumbers.l(maxdeg+3:end,:)=[];
md.solidearth.lovenumbers.l(1:mindeg,:)=0;
md.solidearth.lovenumbers.th=tht(1:maxdeg+1,:);
md.solidearth.lovenumbers.tk=tkt(1:maxdeg+1,:);
md.solidearth.lovenumbers.tl=tlt(1:maxdeg+1,:);
md.solidearth.lovenumbers.pmtf_colinear=pmtf1;
md.solidearth.lovenumbers.pmtf_ortho=pmtf2;
md.solidearth.lovenumbers.timefreq=time;
md.solidearth.rotational.equatorialmoi=8.0131e37;
md.solidearth.rotational.polarmoi=8.0394e37;
md.solidearth.rotational.angularvelocity=7.292115e-5;

%}}}
%mask:  {{{
oceanmask=ones(md.mesh.numberofvertices,1);
icemask=ones(md.mesh.numberofvertices,1);
icemask(pos)=-1;


md.mask.ice_levelset=icemask;
md.mask.ocean_levelset=oceanmask;
% }}}

%time stepping: 
md.timestepping.start_time=0;
md.timestepping.time_step=1000;
md.timestepping.final_time=11000;

time1=md.timestepping.start_time:md.timestepping.time_step:md.timestepping.final_time;
md.masstransport.spcthickness=repmat(md.masstransport.spcthickness, [1 length(time1)]);
%md.masstransport.spcthickness(1:end-1,3:end)=0;
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
md.miscellaneous.name='test2092';

%Solution parameters
md.cluster.np=3;
md.solidearth.settings.reltol=NaN;
md.solidearth.settings.abstol=1e-3;
md.solidearth.settings.sealevelloading=0;
md.solidearth.settings.grdocean=0;
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
md.solidearth.requested_outputs={'SealevelBarystaticIceLoad'};

% max number of iteration reverted back to 10 (i.e., the original default value)
md.solidearth.settings.maxiter=10;

%eustatic + selfattraction + elastic + rotation run:
md.solidearth.settings.selfattraction=1;
md.solidearth.settings.elastic=1;
md.solidearth.settings.viscous=1;
md.solidearth.settings.rotation=1;
md.solidearth.settings.horiz=0;

md=solve(md,'tr');
clear m1 m2 m3
for i=1:length(time1)-1
m1(i)=md.results.TransientSolution(i).SealevelchangePolarMotionX;
m2(i)=md.results.TransientSolution(i).SealevelchangePolarMotionY;
m3(i)=md.results.TransientSolution(i).SealevelchangePolarMotionZ;
end

%Fields and tolerances to track changes
field_names={'PolarmotionX','PolarmotionY', 'PolarmotionZ'};
field_tolerances={4e-05,2e-05,3e-05};
field_values={m1,m2,m3};

return
subplot(1,2,1)
plot(time1(2:end-1),m1(2:end)*180/pi *1e6/md.timestepping.time_step, time1(2:end-1),m2(2:end)*180/pi*1e6/md.timestepping.time_step);

subplot(1,2,2)
plot(time1(1:end-1),cumsum(m1)*180/pi, time1(1:end-1),cumsum(m2)*180/pi)

