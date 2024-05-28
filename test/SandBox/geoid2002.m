%Test Name: EarthSlr

%mesh earth:
md=model;
md.mesh=gmshplanet('radius',6.371012*10^3,'resolution',700.); %500 km resolution mesh
md.cluster.np=3
md.verbose=verbose('11111111');

%parameterize slr solution:
%slr loading:  {{{
md.slr.deltathickness=zeros(md.mesh.numberofelements,1);
md.slr.sealevel=zeros(md.mesh.numberofvertices,1);
md.dsl.global_average_thermosteric_sea_level_change=[0;0];
md.dsl.sea_surface_height_change_above_geoid=zeros(md.mesh.numberofvertices+1,1);
md.dsl.sea_water_pressure_change_at_sea_floor=zeros(md.mesh.numberofvertices+1,1);
%antarctica
late=sum(md.mesh.lat(md.mesh.elements),2)/3;
longe=sum(md.mesh.long(md.mesh.elements),2)/3;
pos=find(late <-80);
md.slr.deltathickness(pos)=-100;
%greenland
pos=find(late > 70 &  late < 80 & longe>-60 & longe<-30);
md.slr.deltathickness(pos)=-100;

%elastic loading from love numbers:
nlov=101;
md.slr.love_h = love_numbers('h','CM'); md.slr.love_h(nlov+1:end)=[];
md.slr.love_k = love_numbers('k','CM'); md.slr.love_k(nlov+1:end)=[];
md.slr.love_l = love_numbers('l','CM'); md.slr.love_l(nlov+1:end)=[];

%}}}
%mask:  {{{
mask=gmtmask(md.mesh.lat,md.mesh.long);
icemask=ones(md.mesh.numberofvertices,1);
pos=find(mask==0);  icemask(pos)=-1;
pos=find(sum(mask(md.mesh.elements),2)<3);   icemask(md.mesh.elements(pos,:))=-1;
md.mask.ice_levelset=icemask;
md.mask.ocean_levelset=-icemask;

%make sure that the elements that have loads are fully grounded:
pos=find(md.slr.deltathickness);
md.mask.ocean_levelset(md.mesh.elements(pos,:))=1;

%make sure wherever there is an ice load, that the mask is set to ice:
pos=find(md.slr.deltathickness);
md.mask.ice_levelset(md.mesh.elements(pos,:))=-1;
% }}}

md.slr.ocean_area_scaling=0;

%Geometry for the bed, arbitrary: 
md.geometry.bed=-ones(md.mesh.numberofvertices,1);

%Materials: 
md.materials=materials('hydro');

%New stuff
md.slr.spcthickness = NaN(md.mesh.numberofvertices,1);
md.slr.Ngia = zeros(md.mesh.numberofvertices,1);
md.slr.Ugia = zeros(md.mesh.numberofvertices,1);
md.slr.hydro_rate = zeros(md.mesh.numberofvertices,1);

%Miscellaneous
md.miscellaneous.name='test2002';

%Solution parameters
md.slr.reltol=NaN;
md.slr.abstol=1e-3;
md.slr.geodetic=1;

% max number of iteration reverted back to 10 (i.e., the original default value)
md.slr.maxiter=10;

md.slr.requested_outputs= {'default',...
	'SealevelriseDeltathickness','Sealevel','SealevelRSLRate','SealevelriseCumDeltathickness',...
	'SealevelNEsaRate', 'SealevelUEsaRate', 'SealevelNGiaRate', 'SealevelUGiaRate',...
	'SealevelEustaticMask','SealevelEustaticOceanMask','SealevelUEsa','SealevelNEsa'};

md.slr.rigid=1; md.slr.elastic=1;md.slr.rotation=1;
md=solve(md,'Sealevelrise');
