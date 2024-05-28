%Test Name: Sea-Level-Partitions
testagainst2002=0;

%Data paths: {{{
shppath='../Data/shp/';
gshhsshapefile=[shppath 'GSHHS_c_L1-NightlyRun.shp'];
%}}}

%create sealevel model to hold our information:
sl=sealevelmodel();

%Create basins using boundaries from shapefile:
%some projections we'll rely on:  %{{{
proj4326=epsg2proj(4326);
proj3031=epsg2proj(3031);
%}}}
%HemisphereWest: {{{
sl.addbasin(basin('continent','hemispherewest','name','hemispherewest','proj',laea(0,-90),'boundaries',{... %Peru projection 3587
	boundary('shppath',shppath,'shpfilename','HemisphereSplit','proj',proj4326,'orientation','reverse'),...
	boundary('shppath',shppath,'shpfilename','NorthAntarctica','proj',proj3031),...
	boundary('shppath',shppath,'shpfilename','RonneBrunt','proj',proj3031,'orientation','reverse'),...
	boundary('shppath',shppath,'shpfilename','RonneEastSummit','proj',proj3031),...
	boundary('shppath',shppath,'shpfilename','RonneFront','proj',proj3031,'orientation','reverse'),...
	boundary('shppath',shppath,'shpfilename','RonneWestSummit','proj',proj3031),...
	boundary('shppath',shppath,'shpfilename','WestAntarctica2','proj',proj3031,'orientation','reverse'),...
	boundary('shppath',shppath,'shpfilename','SouthAntarctica','proj',proj3031)...
	}));
%}}}
%Ross: {{{
sl.addbasin(basin('continent','antarctica','name','ross','proj',proj3031,'boundaries',{...
	boundary('shppath',shppath,'shpfilename','SouthAntarctica','proj',proj3031),...
	boundary('shppath',shppath,'shpfilename','RossIceShelf','proj',proj3031),...
	boundary('shppath',shppath,'shpfilename','RossWestFront','proj',proj3031),...
	boundary('shppath',shppath,'shpfilename','RossFront','proj',proj3031,'orientation','reverse')...
	}));
%}}}
%HemisphereEast: {{{
sl.addbasin(basin('continent','hemisphereeast','name','hemisphereeast','proj',laea(0,+90),'boundaries',{... %Australian projection lat,long
	boundary('shppath',shppath,'shpfilename','HemisphereSplit','proj',proj4326),...
	boundary('shppath',shppath,'shpfilename','SouthAntarctica','proj',proj3031),...
	boundary('shppath',shppath,'shpfilename','RossFront','proj',proj3031),...
	boundary('shppath',shppath,'shpfilename','RossWestFront','proj',proj3031),...
	boundary('shppath',shppath,'shpfilename','EastAntarctica2','proj',proj3031,'orientation','reverse'),...
	boundary('shppath',shppath,'shpfilename','NorthAntarctica','proj',proj3031)...
	}));
%}}}
%Antarctica excluding Ronne: {{{
sl.addbasin(basin('continent','antarctica','name','antarctica-grounded','proj',proj3031,'boundaries',{...
	boundary('shppath',shppath,'shpfilename','NorthAntarctica','proj',proj3031),...
	boundary('shppath',shppath,'shpfilename','EastAntarctica2','proj',proj3031),...
	boundary('shppath',shppath,'shpfilename','RossWestFront','proj',proj3031),...
	boundary('shppath',shppath,'shpfilename','RossIceShelf','proj',proj3031,'orientation','reverse'),...
	boundary('shppath',shppath,'shpfilename','SouthAntarctica','proj',proj3031),...
	boundary('shppath',shppath,'shpfilename','WestAntarctica2','proj',proj3031),...
	boundary('shppath',shppath,'shpfilename','RonneWestSummit','proj',proj3031),...
	boundary('shppath',shppath,'shpfilename','RonneIceShelf','proj',proj3031),...
	boundary('shppath',shppath,'shpfilename','RonneEastSummit','proj',proj3031),...
	boundary('shppath',shppath,'shpfilename','RonneBrunt','proj',proj3031)...
	}));
%}}}
%Ronne: {{{
sl.addbasin(basin('continent','antarctica','name','ronne','proj',proj3031,'boundaries',{...
	boundary('shppath',shppath,'shpfilename','RonneWestSummit','proj',proj3031),...
	boundary('shppath',shppath,'shpfilename','RonneIceShelf','proj',proj3031),...
	boundary('shppath',shppath,'shpfilename','RonneEastSummit','proj',proj3031),...
	boundary('shppath',shppath,'shpfilename','RonneFront','proj',proj3031,'orientation','reverse')...
	}));
%}}}

%Meshing
%Go through basins and mesh:  %{{{
%meshing parameters:  {{{
hmin=500; hmax=700; hmin=hmin*1000; hmax=hmax*1000;
tolerance=100; %tolerance of 100m on Earth position when merging 3d meshes
threshold=5;
defaultoptions={'KeepVertices',0,'MaxCornerAngle',0.0000000001,'NoBoundaryRefinement',1};
alreadyloaded=0;
%}}}
for ind=sl.basinindx('basin','all'),

	bas=sl.basins{ind};
	disp(sprintf('Meshing basin %s\n',bas.name));

	%recover basin domain:
	domain=bas.contour();

	%recover coastline inside basin, using GSHHS_c_L1. It's a lat/long file, hence epsg 4326
	coastline=bas.shapefilecrop('shapefile',gshhsshapefile,'epsgshapefile',4326,'threshold',threshold);

	%mesh:
	md=bamg(model,'domain',domain,'subdomains',coastline,'hmin',hmin,'hmax',hmax,defaultoptions{:});
	%plotmodel(md,'data','mesh');pause(1);

	%miscellaneous:
	md.mesh.proj=bas.proj;
	md.miscellaneous.name=bas.name;

	%recover mask where we have land:
	md.private.bamg.landmask=double(md.private.bamg.mesh.Triangles(:,4)>=1);

	%vertex connectivity:
	md.mesh.vertexconnectivity=NodeConnectivity(md.mesh.elements,md.mesh.numberofvertices);

	%add model to sl icecaps:
	sl.addicecap(md);
end
%}}}

%Parameterization:
%Parameterize ice sheets : {{{
for ind=sl.basinindx('continent',{'antarctica'}),
	disp(sprintf('Parameterizing basin %s\n', sl.icecaps{ind}.miscellaneous.name));

	md=sl.icecaps{ind};
	bas=sl.basins{ind};
	%masks :  %{{{
	%ice levelset from domain outlines:
	md.mask.ice_levelset=-ones(md.mesh.numberofvertices,1);

	if bas.isnameany('antarctica-grounded'),
		md.mask.ocean_levelset=ones(md.mesh.numberofvertices,1);
	end
	if bas.isnameany('ronne','ross'),
		md.mask.ocean_levelset=-ones(md.mesh.numberofvertices,1);
	end
	%}}}
	%latlong:  % {{{
	[md.mesh.long,md.mesh.lat]=gdaltransform(md.mesh.x,md.mesh.y,md.mesh.proj,'EPSG:4326');
	%}}}
	%geometry: {{{
	if bas.iscontinentany('antarctica'),
		di=md.materials.rho_ice/md.materials.rho_water;

		disp('      reading bedrock');
		md.geometry.bed=-ones(md.mesh.numberofvertices,1);
		md.geometry.base=md.geometry.bed;
		md.geometry.thickness=1000*ones(md.mesh.numberofvertices,1);
		md.geometry.surface=md.geometry.bed+md.geometry.thickness;

	end % }}}
	%Slc: {{{
	if bas.iscontinentany('antarctica'),
		if testagainst2002,
			% TODO: Check if the following works as expected: 'pos' is empty, so nothing is assigned to 'md.masstransport.spcthickness(pos)'
			md.masstransport.spcthickness=zeros(md.mesh.numberofvertices,1);
			%antarctica
			late=sum(md.mesh.lat(md.mesh.elements),2)/3;
			longe=sum(md.mesh.long(md.mesh.elements),2)/3;
			pos=find(late <-85);
			ratio=0.225314032985172/0.193045366574523;
			%ratio=   1.276564103522540/.869956;
			md.masstransport.spcthickness(md.mesh.elements(pos,:))= md.masstransport.spcthickness(md.mesh.elements(pos,:))-100*ratio;
		else
			in_fileID=fopen('../Data/AIS_delH_trend.txt', 'r');
			delH=textscan(in_fileID, '%f %f %f');
			fclose(in_fileID);
			longAIS=delH{:,1};
			latAIS=delH{:,2};
			delHAIS=delH{:,3};
			% points=[longAIS,latAIS];
			% index=delaunayn(points);
			index=BamgTriangulate(longAIS, latAIS);
			lat=md.mesh.lat;
			long=md.mesh.long+360;
			pos=find(long>360);
			long(pos)=long(pos)-360;
			delHAIS=InterpFromMesh2d(index,longAIS,latAIS,delHAIS,long,lat);
			northpole=find_point(md.mesh.long,md.mesh.lat,0,90);
			delHAIS(northpole)=0;
			md.masstransport.spcthickness=delHAIS/100;
		end

		md.initialization.sealevel=zeros(md.mesh.numberofvertices,1);

		md.dsl.global_average_thermosteric_sea_level=[0;0];
		md.dsl.sea_surface_height_above_geoid=zeros(md.mesh.numberofvertices+1,1);
		md.dsl.sea_water_pressure_at_sea_floor=zeros(md.mesh.numberofvertices+1,1);

	end %}}}
	% material properties: {{{
	md.materials=materials('hydro');
	%}}}
	%diverse: {{{
	md.miscellaneous.name=bas.name;
	% }}}

	sl.icecaps{ind}=md;
end
%}}}

% ParameterizeContinents {{{
for ind=sl.basinindx('continent',{'hemisphereeast','hemispherewest'}),
	disp(sprintf('Masks for basin %s\n', sl.icecaps{ind}.miscellaneous.name));
	md=sl.icecaps{ind};
	bas=sl.basins{ind};

	%recover lat,long:
	[md.mesh.long,md.mesh.lat]=gdaltransform(md.mesh.x,md.mesh.y,md.mesh.proj,'EPSG:4326');

	%mask:  %{{{
	%Figure out mask from initial mesh: deal with land and ocean masks (land
	%includes grounded ice).  %{{{
	%first, transform land element mask into vertex-driven one
	land=md.private.bamg.landmask;
	land_mask=-ones(md.mesh.numberofvertices,1);

	landels=find(land);
	land_mask(md.mesh.elements(landels,:))=1;

	% Gothrough edges of each land element
	connectedels=md.mesh.elementconnectivity(landels,:);
	connectedisonocean=~land(connectedels);
	sumconnectedisonocean=sum(connectedisonocean,2);

	%figure out which land elements are connected to the ocean:
	landelsconocean=landels(find(sumconnectedisonocean));

	ind1=[md.mesh.elements(landelsconocean,1);
	md.mesh.elements(landelsconocean,2);
	md.mesh.elements(landelsconocean,3)];
	ind2=[md.mesh.elements(landelsconocean,2);
	md.mesh.elements(landelsconocean,3);
	md.mesh.elements(landelsconocean,1)];

	%edge ind1 and ind2:
	for i=1:length(ind1),
		els1=md.mesh.vertexconnectivity(ind1(i),1: md.mesh.vertexconnectivity(ind1(i),end));
		els2=md.mesh.vertexconnectivity(ind2(i),1: md.mesh.vertexconnectivity(ind2(i),end));
		els=intersect(els1,els2);

		if length(find(land(els)))==1,
			%this edge is on the beach, 0 the edge:
			land_mask(ind1(i))=0;
			land_mask(ind2(i))=0;
		end
	end

	md.mask.ocean_levelset=land_mask;
	md.mask.ice_levelset=ones(md.mesh.numberofvertices,1);   %if there are glaciers, we'll adjust

	if testagainst2002,
		% {{{
		%greenland
		pos=find(md.mesh.lat > 70 & md.mesh.lat < 80 & md.mesh.long>-60 & md.mesh.long<-30);
		md.mask.ice_levelset(pos)=-1;
		% }}}
	end
	% }}}
	%}}}
	%slc loading/calibration:  {{{
	md.masstransport.spcthickness=zeros(md.mesh.numberofvertices,1);

	if testagainst2002,
		% {{{
		%greenland
		late=sum(md.mesh.lat(md.mesh.elements),2)/3;
		longe=sum(md.mesh.long(md.mesh.elements),2)/3;
		pos=find(late > 70 &  late < 80 & longe>-60 & longe<-30);
		ratio=.3823/.262344;
		md.masstransport.spcthickness(md.mesh.elements(pos,:))= md.masstransport.spcthickness(md.mesh.elements(pos,:))-100*ratio;
		%md.masstransport.spcthickness(pos)=-100*ratio;

		%correct mask:
		md.mask.ice_levelset(md.mesh.elements(pos,:))=-1;
		% }}}
	else
		delH=textread('../Data/GIS_delH_trend.txt');
		longGIS=delH(:,1);
		latGIS=delH(:,2);
		delHGIS=delH(:,3);
		index=BamgTriangulate(longGIS, latGIS);
		lat=md.mesh.lat;
		long=md.mesh.long+360;
		pos=find(long>360);
		long(pos)=long(pos)-360;
		delHGIS=InterpFromMeshToMesh2d(index,longGIS,latGIS,delHGIS,long,lat);

		delH=textread('../Data/GLA_delH_trend_15regions.txt');
		longGLA=delH(:,1);
		latGLA=delH(:,2);
		delHGLA=sum(delH(:,3:end),2);
		index=BamgTriangulate(longGLA, latGLA);
		lat=md.mesh.lat;
		long=md.mesh.long+360;
		pos=find(long>360);
		long(pos)=long(pos)-360;
		delHGLA=InterpFromMeshToMesh2d(index,longGLA,latGLA,delHGLA,long,lat);

		pos=find(delHGIS);
		md.masstransport.spcthickness(pos) = md.masstransport.spcthickness(pos)-delHGIS(pos)/100;
		pos=find(delHGLA);
		md.masstransport.spcthickness(pos)= md.masstransport.spcthickness(pos)-delHGLA(pos)/100;

		%adjust mask accordingly:
		pos=find(md.masstransport.spcthickness);
		md.mask.ice_levelset(pos)=-1;
		md.mask.ocean_levelset(pos)=1;
	end

	md.initialization.sealevel=zeros(md.mesh.numberofvertices,1);

	md.dsl.global_average_thermosteric_sea_level=[0;0];
	%md.dsl.steric_rate=(1.1+.38)*ones(md.mesh.numberofvertices,1); %steric + water storage.
	md.dsl.sea_surface_height_above_geoid=zeros(md.mesh.numberofvertices+1,1);
	md.dsl.sea_water_pressure_at_sea_floor=zeros(md.mesh.numberofvertices+1,1);

	%}}}
	%geometry:  {{{
	di=md.materials.rho_ice/md.materials.rho_water;
	md.geometry.bed=-ones(md.mesh.numberofvertices,1);
	md.geometry.base=md.geometry.bed;
	md.geometry.thickness=1000*ones(md.mesh.numberofvertices,1);
	md.geometry.surface=md.geometry.bed+md.geometry.thickness;
	% }}}
	%materials:  {{{
	md.materials=materials('hydro');
	% }}}
	sl.icecaps{ind}=md;
end
% }}}

%%Assemble Earth in 3D {{{

%parameters:
plotting=0;
tolerance=100;
loneedgesdetect=0;

%create Earth model by concatenating all the icecaps in 3D:
sl.caticecaps('tolerance',tolerance,'loneedgesdetect',loneedgesdetect);

%figure out how each icecap's mesh connects to the larger Earth mesh:
sl.intersections('force',1);

%figure out connectivity:
disp('Mesh connectivity');
sl.earth.mesh.vertexconnectivity=NodeConnectivity(sl.earth.mesh.elements,sl.earth.mesh.numberofvertices);

%areas:
disp('Mesh nodal areas');
sl.earth.mesh.area=averaging(sl.earth,GetAreas3DTria(sl.earth.mesh.elements,sl.earth.mesh.x,sl.earth.mesh.y,sl.earth.mesh.z),4);

%transfer a list of fields from each icecap and continent back to Earth:
sl.transfer('mask.ice_levelset');
sl.transfer('mask.ocean_levelset');
sl.transfer('geometry.bed');
sl.transfer('geometry.surface');
sl.transfer('geometry.thickness');
sl.transfer('geometry.base');
sl.transfer('mesh.lat');
sl.transfer('mesh.long');
sl.transfer('masstransport.spcthickness');
sl.transfer('initialization.sealevel');
sl.transfer('dsl.sea_surface_height_above_geoid');
sl.transfer('dsl.sea_water_pressure_at_sea_floor');

%radius:
sl.earth.mesh.r=sqrt(sl.earth.mesh.x.^2+sl.earth.mesh.y.^2+sl.earth.mesh.z.^2);

%check on the mesh transitions: {{{
plotting=0;
if plotting,
	flags=ones(sl.earth.mesh.numberofelements,1);
	for i=1:length(sl.eltransitions)
		flags(sl.eltransitions{i})=i;
	end
	plotmodel(sl.earth,'data',flags,'shading','faceted','coastline','on','coast_color','g')
end
%}}}}

% }}}
%Solve Sea-level eqEricuation on Earth only:  {{{
md=sl.earth; %we don't do computations on ice sheets or land.

%Materials:
md.materials=materials('hydro');

%elastic loading from love numbers:
md.solidearth.lovenumbers=lovenumbers('maxdeg',100);
md.solidearth.settings.ocean_area_scaling = 0;

%Miscellaneous
md.miscellaneous.name='test2004';

%New stuff
md.dsl.global_average_thermosteric_sea_level=[1.1+.38;0]; %steric + water storage AR5.

%Solution parameters
md.solidearth.settings.reltol=NaN;
md.solidearth.settings.abstol=1e-3;
md.solidearth.settings.sealevelloading=1;
md.solidearth.settings.isgrd=1;
md.solidearth.settings.ocean_area_scaling=0;
md.solidearth.settings.grdmodel=1;
md.timestepping.time_step=1;

%Physics:
md.transient.issmb=0;
md.transient.isstressbalance=0;
md.transient.isthermal=0;
md.transient.ismasstransport=1;
md.transient.isslc=1;

%Initializations:
md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
md.basalforcings.floatingice_melting_rate=zeros(md.mesh.numberofvertices,1);
md.initialization.vx=zeros(md.mesh.numberofvertices,1);
md.initialization.vy=zeros(md.mesh.numberofvertices,1);
md.initialization.sealevel=zeros(md.mesh.numberofvertices,1);
md.initialization.bottompressure=zeros(md.mesh.numberofvertices,1);
md.initialization.dsl=zeros(md.mesh.numberofvertices,1);
md.initialization.str=0;
md.smb.mass_balance=zeros(md.mesh.numberofvertices,1);

%max number of iterations reverted back to 10 (i.e. the original default value)
md.solidearth.settings.maxiter=10;

%eustatic run:
md.solidearth.settings.selfattraction=0;
md.solidearth.settings.elastic=0;
md.solidearth.settings.rotation=0;
md.solidearth.settings.viscous=0;
md.solidearth.requested_outputs= {'default',...
	'DeltaIceThickness','Sealevel','Bed',...
	'SealevelBarystaticIceMask','SealevelBarystaticOceanMask'};
md=solve(md,'Transient');
Seustatic=md.results.TransientSolution.Sealevel;

%eustatic + selfattraction run:
md.solidearth.settings.selfattraction=1;
md.solidearth.settings.elastic=0;
md.solidearth.settings.rotation=0;
md.solidearth.settings.viscous=0;
md=solve(md,'Transient');
Sselfattraction=md.results.TransientSolution.Sealevel;

%eustatic + selfattraction + elastic run:
md.solidearth.settings.selfattraction=1;
md.solidearth.settings.elastic=1;
md.solidearth.settings.rotation=0;
md.solidearth.settings.viscous=0;
md=solve(md,'Transient');
Selastic=md.results.TransientSolution.Sealevel;

%eustatic + selfattraction + elastic + rotation run:
md.solidearth.settings.selfattraction=1;
md.solidearth.settings.elastic=1;
md.solidearth.settings.rotation=1;
md.solidearth.settings.viscous=0;
md=solve(md,'Transient');
Srotation=md.results.TransientSolution.Sealevel;
%}}}

%Fields and tolerances to track changes
field_names     ={'Eustatic','Rigid','Elastic','Rotation'};
field_tolerances={1e-13,1e-13,1e-13,1e-13};
field_values={Seustatic,Sselfattraction,Selastic,Srotation};
