%Test Name: External_OfflineSolidearthSolution

%mesh earth:
md=model;
load ../Data/SlcTestMesh.mat;
md.mesh=SlcMesh; %700 km resolution mesh

%Geometry for the bed, arbitrary
md.geometry.bed=-ones(md.mesh.numberofvertices,1);
md.geometry.base=md.geometry.bed;
md.geometry.thickness=zeros(md.mesh.numberofvertices,1);
md.geometry.surface=md.geometry.bed+md.geometry.thickness;

%parameterize solidearth solution:
late=md.mesh.lat;
longe=md.mesh.long;
time=0:0.5:5;
%The offline solution pattern is a degree (2,1) spherical harmonic
md.solidearth.external=offlinesolidearthsolution;
md.solidearth.external.displacementup=.5*sind(late).*cosd(late).*cosd(longe) .*time;
md.solidearth.external.displacementup(end+1,:)=time;
md.solidearth.external.geoid=-.1*sind(late).*cosd(late).*sind(longe) .*time;
md.solidearth.external.geoid(end+1,:)=time;
md.solidearth.external.displacementeast=late .*time;
md.solidearth.external.displacementeast(end+1,:)=time;
md.solidearth.external.displacementnorth=longe .*time;
md.solidearth.external.displacementnorth(end+1,:)=time;
md.initialization.sealevel=zeros(md.mesh.numberofvertices,1);


%mask:  {{{
mask=gmtmask(md.mesh.lat,md.mesh.long);
md.mask.ice_levelset=ones(md.mesh.numberofvertices,1);
md.mask.ocean_levelset=ones(md.mesh.numberofvertices,1);
% }}}

%time stepping: 
md.timestepping.start_time=time(1);
md.timestepping.time_step=time(2)-time(1);
md.timestepping.final_time=time(end);

%Materials: 
md.materials=materials('hydro');

%Miscellaneous
md.miscellaneous.name='test2007';

%Solution parameters
md.cluster.np=3;
md.solidearth.settings.isgrd=0;
md.solidearth.settings.horiz=1;

%Physics: bary
md.transient.issmb=0; 
md.transient.isstressbalance=0;
md.transient.isthermal=0;
md.transient.ismasstransport=0;
md.transient.isslc=1;
md.solidearth.requested_outputs={'Sealevel', 'Bed', 'BedEast', 'BedNorth'};

%eustatic run:
md=solve(md,'Transient');

for i=length(time)-1;
Geoid=md.results.TransientSolution(i).Sealevel;
BedUp=md.results.TransientSolution(i).Bed;
BedEast=md.results.TransientSolution(i).BedEast;
BedNorth=md.results.TransientSolution(i).BedNorth;
end

%Fields and tolerances to track changes
field_names={'Geoid','BedUp','BedEast','BedNorth'};
field_tolerances={1e-13,1e-13,1e-13,1e-13};
field_values={Geoid,BedUp,BedEast,BedNorth};
