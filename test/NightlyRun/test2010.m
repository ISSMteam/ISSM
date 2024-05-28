%Test Name: MomentOfIntertia

%mesh earth:
md=model;
load ../Data/SlcTestMesh.mat;
md.mesh=SlcMesh; %700 km resolution mesh

%Geometry for the bed, arbitrary thickness of 100:
md.geometry.bed=-ones(md.mesh.numberofvertices,1);
md.geometry.base=md.geometry.bed;
md.geometry.thickness=100*ones(md.mesh.numberofvertices,1);
md.geometry.surface=md.geometry.bed+md.geometry.thickness;

%parameterize slc solution:
%solidearth loading:  {{{
md.masstransport.spcthickness=[md.geometry.thickness;0];
md.smb.mass_balance=zeros(md.mesh.numberofvertices,1);

xe=md.mesh.x(md.mesh.elements)*[1;1;1]/3;
ye=md.mesh.y(md.mesh.elements)*[1;1;1]/3;
ze=md.mesh.z(md.mesh.elements)*[1;1;1]/3;
re=sqrt(xe.^2+ye.^2+ze.^2);

late=asind(ze./re);
longe=atan2d(ye,xe);
%greenland
pos=find(late>60 & late<90 & longe>-75 & longe<-15);
md.masstransport.spcthickness(md.mesh.elements(pos,:))= md.masstransport.spcthickness(md.mesh.elements(pos,:))-100;
posice=pos;

%elastic loading from love numbers:
md.solidearth.lovenumbers=lovenumbers('maxdeg',100);

%}}}
%mask:  {{{
mask=gmtmask(md.mesh.lat,md.mesh.long);
icemask=ones(md.mesh.numberofvertices,1);
icemask(md.mesh.elements(posice,:))=-0.5;

oceanmask=-ones(md.mesh.numberofvertices,1);
pos=find(mask==0); oceanmask(pos)=1; icemask(~pos)=1;

md.mask.ice_levelset=icemask;
md.mask.ocean_levelset=oceanmask;

% use model representation of ocean area (not the true area)
md.solidearth.settings.ocean_area_scaling = 0;

%materials
md.initialization.temperature=273.25*ones(md.mesh.numberofvertices,1);
md.initialization.sealevel=zeros(md.mesh.numberofvertices,1);
md.initialization.str=0;

md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
md.basalforcings.floatingice_melting_rate=zeros(md.mesh.numberofvertices,1);
md.initialization.vx=zeros(md.mesh.numberofvertices,1);
md.initialization.vy=zeros(md.mesh.numberofvertices,1);

%Miscellaneous
md.miscellaneous.name='test2010';

%Solution parameters
md.solidearth.settings.reltol=NaN;
md.solidearth.settings.abstol=1e-3;
md.solidearth.settings.sealevelloading=0;
md.solidearth.settings.grdocean=1;
md.solidearth.settings.isgrd=1;
md.solidearth.settings.ocean_area_scaling=0;
md.solidearth.settings.grdmodel=1;
md.solidearth.settings.horiz=1;
md.solidearth.requested_outputs={'Sealevel','SealevelBarystaticIceArea','SealevelBarystaticIceLoad','SealevelBarystaticIceMask','SealevelBarystaticIceLatbar','SealevelBarystaticIceLongbar'};

%Physics: 
md.transient.issmb=0; 
md.transient.isstressbalance=0;
md.transient.isthermal=0;
md.transient.ismasstransport=1;
md.transient.isslc=1;

md.timestepping.start_time=0;
md.timestepping.time_step=1;
md.timestepping.final_time=1;

%eustatic + selfattraction + elastic + rotation run:
md.solidearth.settings.selfattraction=1;
md.solidearth.settings.elastic=1;
md.solidearth.settings.rotation=1;
md.solidearth.settings.viscous=0;
md.cluster=generic('name',oshostname(),'np',3);
%md.verbose=verbose('111111111');
md=solve(md,'Transient');

moi_p=md.solidearth.rotational.polarmoi;
moi_e=md.solidearth.rotational.equatorialmoi;
tide_love_k2=md.solidearth.lovenumbers.tk(3);
load_love_k2=md.solidearth.lovenumbers.k(3);
tide_love_k2secular=md.solidearth.lovenumbers.tk2secular;
% uncomment following 2 lines for
eus=md.results.TransientSolution.Bslc;
slc=md.results.TransientSolution.Sealevel;
moixz=md.results.TransientSolution.SealevelchangePolarMotionX / (1/(1-tide_love_k2/tide_love_k2secular) * (1+load_love_k2)/(moi_p-moi_e) );
moiyz=md.results.TransientSolution.SealevelchangePolarMotionY / (1/(1-tide_love_k2/tide_love_k2secular) * (1+load_love_k2)/(moi_p-moi_e) );
moizz=md.results.TransientSolution.SealevelchangePolarMotionZ / ( -(1+load_love_k2)/moi_p);

areaice=md.results.TransientSolution.SealevelBarystaticIceArea;
areaice(isnan(areaice))=0;
loadice=md.results.TransientSolution.SealevelBarystaticIceLoad;
rad_e = md.solidearth.planetradius;

lat=md.results.TransientSolution.SealevelBarystaticIceLatbar*pi/180;
lon=md.results.TransientSolution.SealevelBarystaticIceLongbar*pi/180;
moi_xz = sum(-loadice.*rad_e^2.*sin(lat).*cos(lat).*cos(lon));
moi_yz = sum(-loadice.*rad_e^2.*sin(lat).*cos(lat).*sin(lon));
moi_zz = sum(-loadice.*rad_e^2.*(-1.0/3.0+sin(lat).^2));
theoretical_value_check=[moixz/moi_xz moiyz/moi_yz moizz/moi_zz]; % should yield [1.0 1.0 1.0]
% }}}

%Fields and tolerances to track changes
field_names     ={'eus','slc','moixz','moiyz','moizz'};
field_tolerances={1e-13,1e-13,1e-13,1e-13,1e-13};
field_values={eus,slc,moixz,moiyz,moizz};

