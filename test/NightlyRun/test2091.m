%Test Name: ViscoElasticEarthSlc

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
rad=re*pi*[linspace(1e-3,1,3)/18];
for i=1:length(rad);
co(i).x=cos(phi2')*rad(i);
co(i).y=sin(phi2')*rad(i);co(i).y(end)=co(i).y(1);
co(i).nods=length(co(i).x);co(i).closed=1;co(i).density=1;co(i).name='Icedisklimit';co(i).Geometry='Polygon';
co(i).BoundingBox=[min(co(i).x) min(co(i).y); max(co(i).x) max(co(i).y)];
end

profile_colat=(0:1:179)';
profile_x=profile_colat/180*pi*re;
profile_y=profile_x*0;
requiredprofile=true(length(profile_colat),1);
geometry_x=dom.x;
geometry_y=dom.y;
for i=1:length(co)
	geometry_x=[geometry_x;co(i).x];
	geometry_y=[geometry_y;co(i).y];
end
for i=1:length(profile_colat)
	if min((geometry_x-profile_x(i)).^2+(geometry_y-profile_y(i)).^2)<1.
		requiredprofile(i)=false;
	end
end
requiredvertices=[profile_x(requiredprofile) profile_y(requiredprofile)];

defaultoptions={'KeepVertices',0,'MaxCornerAngle',0.0000000001,'NoBoundaryRefinement',1,'RequiredVertices',requiredvertices};
md2d=bamg(model,'domain',dom,'subdomains',co,'hmin',100e3,'hmax',10000e3,defaultoptions{:});

profile_indices=zeros(length(profile_colat),1);
for i=1:length(profile_colat)
	dist=(md2d.mesh.x-profile_x(i)).^2+(md2d.mesh.y-profile_y(i)).^2;
	[distmin,profile_indices(i)]=min(dist);
	if(distmin>1.) 
        error('Could not recover required profile vertex in test2091');
    end
end

colat=sqrt(md2d.mesh.x.^2+md2d.mesh.y.^2)/re;
lon=atan2(md2d.mesh.y,md2d.mesh.x);

x=sin(colat).*cos(lon)*re;
y=sin(colat).*sin(lon)*re;
z=cos(colat)*re;

md.mesh.x=x;md.mesh.y=y;md.mesh.z=z;md.mesh.lat=90-colat*180/pi;md.mesh.long=lon*180/pi;
md.mesh.r=sqrt(x.^2+y.^2+z.^2);
md.mesh.numberofelements=md2d.mesh.numberofelements;md.mesh.numberofvertices=md2d.mesh.numberofvertices;
md.mesh.elements=md2d.mesh.elements;
md.mesh.vertexconnectivity=NodeConnectivity(md.mesh.elements,md.mesh.numberofvertices);
md.mesh.area=GetAreas3DTria(md.mesh.elements,md.mesh.x,md.mesh.y,md.mesh.z);


%md.mesh=gmshplanet('radius',6.371012*10^3,'resolution',500.);

%Geometry for the bed, arbitrary thickness of 1000: 
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
%arctic
pos=find(late>=80);
pos1=find(md.mesh.lat>=80);
md.masstransport.spcthickness(:)=0;
md.masstransport.spcthickness(pos1)=md.masstransport.spcthickness(pos1) +1500*sqrt((cosd(90-md.mesh.lat(pos1))-cosd(10))./(1-cosd(10)));
%md.masstransport.spcthickness(pos1)=md.masstransport.spcthickness(pos1) +1000;
posgre=pos;

%visco-elastic loading from love numbers that we load ourselves. We 
%still use the lovenumbers constructor to initialize the fields: 
load ../Data/lnb_temporal.mat
%md.solidearth.lovenumbers=lovenumbers('maxdeg',1000);
%md.solidearth.lovenumbers.timefreq=[0];

maxdeg=129;
mindeg=1;

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
%md.solidearth.lovenumbers.h=repmat(md.solidearth.lovenumbers.h,1,100);
%md.solidearth.lovenumbers.k=repmat(md.solidearth.lovenumbers.k,1,100);
%md.solidearth.lovenumbers.l=repmat(md.solidearth.lovenumbers.l,1,100);
md.solidearth.lovenumbers.th=tht(1:maxdeg+1,:);
md.solidearth.lovenumbers.tk=tkt(1:maxdeg+1,:);
md.solidearth.lovenumbers.tl=tlt(1:maxdeg+1,:);
md.solidearth.lovenumbers.pmtf_colinear=pmtf1;
md.solidearth.lovenumbers.pmtf_ortho=pmtf2;
md.solidearth.lovenumbers.timefreq=time;

%}}}
%mask:  {{{
%mask=gmtmask(md.mesh.lat,md.mesh.long);
oceanmask=ones(md.mesh.numberofvertices,1);

icemask=ones(md.mesh.numberofvertices,1);
icemask(pos1)=-1;


md.mask.ice_levelset=icemask*1e-3;
md.mask.ocean_levelset=oceanmask;
% }}}

%time stepping: 
md.timestepping.start_time=0;
md.timestepping.time_step=1000;
md.timestepping.final_time=12000;

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
md.miscellaneous.name='test2091';

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
md.solidearth.settings.degacc=.01;
md.solidearth.settings.viscoussampling=20;

%Physics:
md.transient.issmb=0; 
md.transient.isstressbalance=0;
md.transient.isthermal=0;
md.transient.ismasstransport=1;
md.transient.isslc=1;
md.solidearth.requested_outputs={'SealevelGRD','BedGRD', 'BedNorthGRD', 'BedEastGRD', 'SealevelBarystaticIceLoad', 'SealevelBarystaticIceWeights', 'SealevelBarystaticIceArea', 'SealevelBarystaticIceMask', 'SealevelBarystaticIceLatbar', 'SealevelBarystaticIceLongbar'};
md.settings.results_on_nodes = {'SealevelBarystaticIceWeights'};

% max number of iteration reverted back to 10 (i.e., the original default value)
md.solidearth.settings.maxiter=10;

%eustatic + selfattraction + elastic + rotation run:
md.solidearth.settings.selfattraction=1;
md.solidearth.settings.elastic=1;
md.solidearth.settings.viscous=1;
md.solidearth.settings.rotation=0;
md.solidearth.settings.horiz=1;

md=solve(md,'tr');


%%validation against spada curves
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
E=cumsum(E,2); %E should theoretically be 0, in practice there is a residual due to the resolution and imperfect symmetry of the ice load. This can serve as an accuracy check of the horizontals

%Archive values at required profile vertices rather than raw node fields. The
%BAMG mesh for this circular-cap benchmark can have platform-dependent node counts.
S=S(profile_indices,:);
B=B(profile_indices,:);
H=H(profile_indices,:);

%Fields and tolerances to track changes
field_names={'Sealevel','Bed', 'Horizontals'};
field_tolerances={1e-13,1e-13,1e-13};
field_values={S,B,H};

return
%%
%Spada benchmark curves
[ali,ia,~]=unique(90-md.mesh.lat);
fid=fopen('../Data/SpadaBenchmark/T1-2/GS-disp_t0_d1_cap.dat');
 dat=textscan(fid,'%f','headerlines',16);
 dat2=reshape(dat{1}, [5 180*10+2]);
 ur(:,1)=dat2(2,:);
 ut(:,1)=dat2(3,:);
 g(:,1)=dat2(5,:);
 fid=fopen('../Data/SpadaBenchmark/T1-2/GS-disp_t1_d1_cap.dat');
 dat=textscan(fid,'%f','headerlines',16);
 dat2=reshape(dat{1}, [5 180*10+2]);
 ur(:,end+1)=dat2(2,:);
 ut(:,end+1)=dat2(3,:);
 g(:,end+1)=dat2(5,:);
 fid=fopen('../Data/SpadaBenchmark/T1-2/GS-disp_t2_d1_cap.dat');
 dat=textscan(fid,'%f','headerlines',16);
 dat2=reshape(dat{1}, [5 180*10+2]);     
 ur(:,end+1)=dat2(2,:);
 ut(:,end+1)=dat2(3,:);
 g(:,end+1)=dat2(5,:);
 fid=fopen('../Data/SpadaBenchmark/T1-2/GS-disp_t5_d1_cap.dat');
 dat=textscan(fid,'%f','headerlines',16);
 dat2=reshape(dat{1}, [5 180*10+2]);     
 ur(:,end+1)=dat2(2,:);
 ut(:,end+1)=dat2(3,:);
 g(:,end+1)=dat2(5,:);
 fid=fopen('../Data/SpadaBenchmark/T1-2/GS-disp_t10_d1_cap.dat');
 dat=textscan(fid,'%f','headerlines',16);
 dat2=reshape(dat{1}, [5 180*10+2]);     
 ur(:,end+1)=dat2(2,:);
 ut(:,end+1)=dat2(3,:);
 g(:,end+1)=dat2(5,:);

%correction factor for total mass in Spada vs us
fact=3.607171340900778E+018/sum(md.results.TransientSolution(1).SealevelBarystaticIceLoad.*md.results.TransientSolution(1).SealevelBarystaticIceArea);
[ali,ia,~]=unique(90-md.mesh.lat);

clear uri;
clear uti;
clear gi;
for i=1:5
	uri(:,i)=interp1(dat2(1,1:end-1)',ur(1:end-1,i),ali);
	uti(:,i)=interp1(dat2(1,1:end-1)',ut(1:end-1,i),ali);
	gi(:,i)=interp1(dat2(1,1:end-1)',g(1:end-1,i),ali);
end
indB=[1 2 3 6 11];%curves at 0, 1, 2, 5 and 10kyr


return
subplot(2,3,1)
plot(ali,B(ia,indB)*fact,'.',ali,uri)
title('bedrock up');

subplot(2,3,2)
plot(ali,-H(ia,indB)*fact,'.', ali,uti) %minus sign is because unit vector e_North=-e_theta (colatitude direction)
title('bedrock horiz');

subplot(2,3,3)
plot(ali,(S(ia,indB))*fact,'.',ali,gi)
title('Geoid');

subplot(2,3,4)
plot(ali,B(ia,indB)*fact-uri)
title('Bedrock up error');

subplot(2,3,5)
plot(ali,-H(ia,indB)*fact-uti)
title('Bedrock horiz error');

subplot(2,3,6)
plot(ali,(S(ia,indB))*fact-gi)
title('Geoid error');
