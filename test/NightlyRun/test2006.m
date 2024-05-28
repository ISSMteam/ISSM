%Test Name: EarthSlc Dakota Sampling glaciers.

%mesh earth:
md=model;
load ../Data/SlcTestMesh.mat;
md.mesh=SlcMesh; %700 km resolution mesh

%Geometry for the bed, arbitrary thickness of 100: 
md.geometry.bed=zeros(md.mesh.numberofvertices,1);
md.geometry.base=md.geometry.bed;
md.geometry.thickness=100*ones(md.mesh.numberofvertices,1);
md.geometry.surface=md.geometry.bed+md.geometry.thickness;

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
md.masstransport.spcthickness(md.mesh.elements(pos,:))= md.masstransport.spcthickness(md.mesh.elements(pos,:))-100;
posant=pos;
%greenland
pos=find(late>60 & late<90 & longe>-75 & longe<-15);
md.masstransport.spcthickness(md.mesh.elements(pos,:))= md.masstransport.spcthickness(md.mesh.elements(pos,:))-100;
posgre=pos;

%elastic loading from love numbers:
md.solidearth.lovenumbers=lovenumbers('maxdeg',100);

%}}}
%mask:  {{{
mask=gmtmask(md.mesh.lat,md.mesh.long);
oceanmask=-ones(md.mesh.numberofvertices,1);
pos=find(mask==0); oceanmask(pos)=1;

icemask=ones(md.mesh.numberofvertices,1);
icemask(md.mesh.elements(posant))=-1;
icemask(md.mesh.elements(posgre))=-1;

md.mask.ice_levelset=icemask;
md.mask.ocean_levelset=oceanmask;
% }}}

%time stepping: 
md.timestepping.start_time=0;
md.timestepping.time_step=1;
md.timestepping.final_time=10;

%masstransport:
md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
md.basalforcings.floatingice_melting_rate=zeros(md.mesh.numberofvertices,1);
md.initialization.vx=zeros(md.mesh.numberofvertices,1);
md.initialization.vy=zeros(md.mesh.numberofvertices,1);
md.initialization.sealevel=zeros(md.mesh.numberofvertices,1);
md.initialization.str=0;

%Materials: 
md.materials=materials('hydro');

%Miscellaneous
md.miscellaneous.name='test2006';

%Solution parameters
md.cluster.np=3;
md.solidearth.settings.reltol=NaN;
md.solidearth.settings.abstol=1e-3;
md.solidearth.settings.sealevelloading=1;
md.solidearth.settings.isgrd=1;
md.solidearth.settings.ocean_area_scaling=0;
md.solidearth.settings.grdmodel=1;

md.solidearth.settings.selfattraction=1;
md.solidearth.settings.elastic=1;
md.solidearth.settings.rotation=1;
md.solidearth.settings.viscous=0;

%Physics: 
md.transient.issmb=0; 
md.transient.isstressbalance=0;
md.transient.isthermal=0;
md.transient.ismasstransport=1;
md.transient.isslc=1;
md.solidearth.requested_outputs={'Sealevel'};

dh=md.masstransport.spcthickness;
deltathickness=zeros(md.mesh.numberofvertices+1,10);
for i=0:10,
	deltathickness(1:end-1,i+1)=md.geometry.thickness+dh(1:end-1)*i;
end
deltathickness(end,:)=0:1:10;
md.masstransport.spcthickness=deltathickness;

%hack: 
md.geometry.surface=zeros(md.mesh.numberofvertices,1);
md.geometry.thickness=ones(md.mesh.numberofvertices,1);
md.geometry.base=-ones(md.mesh.numberofvertices,1);
md.geometry.bed=md.geometry.base;


%Uncertainty Quantification
%ice sheets {{{
npart=1; nt=1;
partition=-ones(md.mesh.numberofelements,1);
pos=find(late < -80); partition(pos)=0;
pos=find(late>70 & late<80 & longe>-60 & longe<-30); partition(pos)=0;

%variables: 
qmuvar.surfaceload=normal_uncertain('descriptor','scaled_SurfaceloadIceThicknessChange',...
	'mean',1*ones(npart,nt),...
	'stddev',.1*ones(npart,nt),... %10% standard deviation
	'partition',partition,...
	'transient','on',...
	'nsteps',nt);
%}}}

%correlation:
md.qmu.correlation_matrix=[];

%variables final declaration:
md.qmu.variables=struct();
md.qmu.variables.surfaceload=qmuvar.surfaceload;

locations=[1 5 10 15 20];
%responses  % {{{
md.qmu.responses.sealevel1=response_function('descriptor','Outputdefinition1');
md.qmu.responses.sealevel2=response_function('descriptor','Outputdefinition2');
md.qmu.responses.sealevel3=response_function('descriptor','Outputdefinition3');
md.qmu.responses.sealevel4=response_function('descriptor','Outputdefinition4');
md.qmu.responses.sealevel5=response_function('descriptor','Outputdefinition5');

%output definitions: 
for i=1:length(locations),
	ind=locations(i);
	if i==1,
		md.outputdefinition.definitions={nodalvalue('name','SNode','definitionstring','Outputdefinition1', ...
			'model_string','Sealevel','node',ind)}; 
	else
		md.outputdefinition.definitions{end+1}=nodalvalue('name','SNode','definitionstring',['Outputdefinition' num2str(i)], ...
			'model_string','Sealevel','node',ind); 
	end
end 
% }}}

%algorithm:  % {{{
md.qmu.method     =dakota_method('nond_samp');
md.qmu.method(end)=dmeth_params_set(md.qmu.method(end),...
'seed',1234,...
'samples',10,...
'sample_type','random');
md.qmu.output=1; 
%}}}
%parameters % {{{
md.qmu.params.direct=true;
md.qmu.params.interval_type='forward';
md.qmu.params.analysis_driver='matlab';
md.qmu.params.evaluation_scheduling='master';
md.qmu.params.processors_per_evaluation=2;
md.qmu.params.tabular_graphics_data=true; 
md.qmu.isdakota=1;
md.verbose=verbose(0); md.verbose.qmu=1;
% }}}
%qmu statistics %{{{
md.qmu.statistics.nfiles_per_directory=2;
md.qmu.statistics.ndirectories=5;

md.qmu.statistics.method(1).name='Histogram';
md.qmu.statistics.method(1).fields={'Sealevel','BslcIce'};
md.qmu.statistics.method(1).steps=[1:10];
md.qmu.statistics.method(1).nbins=20;

md.qmu.statistics.method(2).name='MeanVariance';
md.qmu.statistics.method(2).fields={'Sealevel','BslcIce'};
md.qmu.statistics.method(2).steps=[1:10];

md.qmu.statistics.method(3).name='SampleSeries';
md.qmu.statistics.method(3).fields={'Sealevel','BslcIce'};
md.qmu.statistics.method(3).steps=[1:10];
md.qmu.statistics.method(3).indices=locations;
%}}}

%run transient dakota solution: 
mds=solve(md,'Transient');

%run without statistics computations:
md.qmu.statistics.method(1).name='None';
md=solve(md,'Transient');

%compare statistics with our own here:
svalues=mds.results.StatisticsSolution(end).SealevelSamples; %all values at locations.

dvalues=zeros(md.qmu.method.params.samples,length(locations));
for i=1:md.qmu.method.params.samples,
	dvalues(i,:)=md.results.dakota.modelresults{i}.TransientSolution(end).Sealevel(locations);
end

samplesnorm=norm(dvalues-svalues,'fro');

%Fields and tolerances to track changes
field_names={'Samples Norm'};
field_tolerances={1e-13};
field_values={samplesnorm};
