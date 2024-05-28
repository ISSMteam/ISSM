%Test Name: EarthSlr

%mesh earth:
md=model;
md.mesh=gmshplanet('radius',6.371012*10^3,'resolution',1000.); %500 km resolution mesh

%parameterize slr solution:
%slr loading:  {{{
md.solidearth.surfaceload.icethicknesschange=zeros(md.mesh.numberofelements,1);
md.solidearth.surfaceload.waterheightchange=zeros(md.mesh.numberofelements,1);
md.solidearth.sealevel=zeros(md.mesh.numberofvertices,1);
md.dsl.global_average_thermosteric_sea_level_change=[0;0];
md.dsl.sea_surface_height_change_above_geoid=zeros(md.mesh.numberofvertices+1,1);
md.dsl.sea_water_pressure_change_at_sea_floor=zeros(md.mesh.numberofvertices+1,1);
%US
late=sum(md.mesh.lat(md.mesh.elements),2)/3;
longe=sum(md.mesh.long(md.mesh.elements),2)/3;

pos=find(late < (37+5) &  late > (37-5) & longe>(-95-10) & longe<(-95+10));
md.solidearth.surfaceload.waterheightchange(pos)=-100;

%elastic loading from love numbers:
md.solidearth.lovenumbers=lovenumbers('maxdeg',100);

%}}}
%mask:  {{{
mask=gmtmask(md.mesh.lat,md.mesh.long);
icemask=ones(md.mesh.numberofvertices,1);
pos=find(mask==0);  icemask(pos)=-1;
pos=find(sum(mask(md.mesh.elements),2)<3);   icemask(md.mesh.elements(pos,:))=-1;
md.mask.ice_levelset=icemask;
md.mask.ocean_levelset=-icemask;

%make sure that the elements that have loads are fully grounded:
pos=find(md.solidearth.surfaceload.waterheightchange);
md.mask.ocean_levelset(md.mesh.elements(pos,:))=1;
md.mask.ice_levelset(md.mesh.elements(pos,:))=1;

% }}}

md.solidearth.settings.ocean_area_scaling=0;

%Geometry for the bed, arbitrary: 
md.geometry.bed=-ones(md.mesh.numberofvertices,1);

%Materials: 
md.materials=materials('hydro');

%Hydro load: 


%Miscellaneous
md.miscellaneous.name='test2002';

%Solution parameters
md.solidearth.settings.reltol=NaN;
md.solidearth.settings.abstol=1e-3;
md.solidearth.settings.computesealevelchange=1;

% max number of iteration reverted back to 10 (i.e., the original default value)
md.solidearth.settings.maxiter=10;

%eustatic run:
md.solidearth.settings.rigid=0; md.solidearth.settings.elastic=0;md.solidearth.settings.rotation=0;
md=solve(md,'Sealevelrise');
Seustatic=md.results.SealevelriseSolution.Sealevel;

%eustatic + rigid run:
md.solidearth.settings.rigid=1; md.solidearth.settings.elastic=0;md.solidearth.settings.rotation=0;
md=solve(md,'Sealevelrise');
Srigid=md.results.SealevelriseSolution.Sealevel;

%eustatic + rigid + elastic run:
md.solidearth.settings.rigid=1; md.solidearth.settings.elastic=1;md.solidearth.settings.rotation=0;
md=solve(md,'Sealevelrise');
Selastic=md.results.SealevelriseSolution.Sealevel;

%eustatic + rigid + elastic + rotation run:
md.solidearth.settings.rigid=1; md.solidearth.settings.elastic=1; md.solidearth.settings.rotation=1;
md=solve(md,'Sealevelrise');
Srotation=md.results.SealevelriseSolution.Sealevel;

%Fields and tolerances to track changes
field_names     ={'Eustatic','Rigid','Elastic','Rotation'};
field_tolerances={1e-13,1e-13,1e-13,1e-13};
field_values={Seustatic,Srigid,Selastic,Srotation};
