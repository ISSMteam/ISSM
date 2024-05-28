%Test Name: EarthSlc_Geometry

step=[1];
if any(step==1)
%mesh earth:
md=model;
load ../Data/SlcTestMesh.mat;
md.mesh=SlcMesh; %700 km resolution mesh

%Geometry for the bed, arbitrary thickness of 1000: 
md.geometry.bed=-ones(md.mesh.numberofvertices,1);
md.geometry.base=md.geometry.bed;
md.geometry.thickness=1000*ones(md.mesh.numberofvertices,1);
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
md.masstransport.spcthickness(md.mesh.elements(:,:))= md.masstransport.spcthickness(md.mesh.elements(:,:))-1000;
posice=pos;

%elastic loading from love numbers:
md.solidearth.lovenumbers=lovenumbers('maxdeg',100);
%}}}
%mask:  {{{
mask=gmtmask(md.mesh.lat,md.mesh.long);
icemask=ones(md.mesh.numberofvertices,1);
icemask(md.mesh.elements(posice,:))=-1;
md.mask.ice_levelset=icemask;
oceanmask=-ones(md.mesh.numberofvertices,1);
pos=find(mask==0); oceanmask(pos)=3;
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
md.miscellaneous.name='test2013';

%Solution parameters
md.solidearth.settings.reltol=NaN;
md.solidearth.settings.abstol=1e-3;
md.solidearth.settings.sealevelloading=0;
md.solidearth.settings.grdocean=0;
md.solidearth.settings.isgrd=1;
md.solidearth.settings.ocean_area_scaling=0;
md.solidearth.settings.grdmodel=1;
md.solidearth.settings.horiz=1;
md.settings.results_on_nodes = {'SealevelBarystaticIceWeights'};
md.solidearth.requested_outputs={'SealevelBarystaticIceLoad', 'SealevelBarystaticIceWeights', 'SealevelBarystaticIceArea', 'SealevelBarystaticIceMask', 'SealevelBarystaticIceLatbar', 'SealevelBarystaticIceLongbar'};

%Physics: 
md.transient.issmb=0;
md.transient.isstressbalance=0;
md.transient.isthermal=0;
md.transient.ismasstransport=1;
md.transient.isslc=1;

md.timestepping.start_time=0;
md.timestepping.time_step=1;
md.timestepping.final_time=1;

%slc:
md.solidearth.settings.selfattraction=1;
md.solidearth.settings.elastic=0;
md.solidearth.settings.rotation=0;
md.solidearth.settings.viscous=0;
md.cluster=generic('name',oshostname(),'np',3);
%md.verbose=verbose('111111111');
md=solve(md,'Transient');

weights=md.results.TransientSolution.SealevelBarystaticIceWeights;
mask=md.results.TransientSolution.SealevelBarystaticIceMask;
loads=md.results.TransientSolution.SealevelBarystaticIceLoad;
areas=md.results.TransientSolution.SealevelBarystaticIceArea;
longbar=md.results.TransientSolution.SealevelBarystaticIceLongbar;
latbar=md.results.TransientSolution.SealevelBarystaticIceLatbar;

loads=loads./areas;

%Fields and tolerances to track changes
field_names     ={'Mask', 'LoadAreas', 'SurfaceLoad', 'LoadWeights','LatitudeLoadBarycenter', 'LongitudeLoadBarycenter'};
field_tolerances={1e-13,1e-13,1e-13,1e-13,1e-13,1e-13};
field_values={mask,areas,loads,weights,latbar, longbar};


end
if any(step==2)
%use this to check visually the load-related inputs

	 plotmodel(md,'data',md.mask.ocean_levelset,'contourlevels', {0},'contouronly',1, 'contourcolor', 'b')  
	 co=get(gca,'children');
	 coo=struct;
	for i=1:length(co)
		coo(i).XData=co(i).XData;
		coo(i).YData=co(i).YData;
		coo(i).ZData=co(i).ZData;
	end	

	 plotmodel(md,'data',md.mask.ice_levelset,'contourlevels', {0},'contouronly',1, 'contourcolor', 'r')  
	 ci=get(gca,'children');
	cii=struct;
	for i=1:length(ci)
		cii(i).XData=ci(i).XData;
		cii(i).YData=ci(i).YData;
		cii(i).ZData=ci(i).ZData;
	end	

	close all

	subplot(2,2,1)
	 for i=1:md.mesh.numberofelements
		patch('XData',md.mesh.x(md.mesh.elements(i,:)), 'YData',md.mesh.y(md.mesh.elements(i,:)), 'ZData',md.mesh.z(md.mesh.elements(i,:)),'facecolor', 'interp', 'facevertexcdata',weights(i,:)')
		hold on
	end

	for i=1:length(coo)
		plot3(coo(i).XData,coo(i).YData,coo(i).ZData, 'color', [0 0 .8],'linewidth',1.5)
	end

	for i=1:length(cii)
		plot3(cii(i).XData,cii(i).YData,cii(i).ZData, 'color', [.7 0 .7],'linewidth',1.5)
	end

	x=md.solidearth.planetradius* (cosd(latbar).*cosd(longbar));
	y=md.solidearth.planetradius* (cosd(latbar).*sind(longbar));
	z=md.solidearth.planetradius* (sind(latbar));
	ind=find(mask>0);
	plot3(x(ind),y(ind),z(ind), 'kx')
	axis tight; title('Load weights & barycenters');set(gca,'Fontsize', 14);

	 plotmodel(md,'data', loads,'subplot',[2 2 2]);title('Average load [kg.m^-2]')
	 plotmodel(md,'data', mask,'subplot',[2 2 3]); title('Phi')
	 plotmodel(md,'data', areas,'subplot',[2 2 4]); title('Load areas [m^2]') 


end

