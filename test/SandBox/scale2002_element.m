%Test Name: Dakota EarthSlr Scaling.

%mesh earth:
md=model;
md.mesh=gmshplanet('radius',6.371012*10^3,'resolution',700.); %500 km resolution mesh
md.cluster.np=2;

%parameterize slr solution:
nsteps=100;
md.solidearth.surfaceload.icethicknesschange=ones(md.mesh.numberofelements+1,nsteps);
md.solidearth.surfaceload.icethicknesschange(end,:)=2000.5:1:2099.5;
md.solidearth.sealevel=zeros(md.mesh.numberofvertices,1);

%antarctica
late=sum(md.mesh.lat(md.mesh.elements),2)/3;
longe=sum(md.mesh.long(md.mesh.elements),2)/3;

posant=find(late <-80);
%md.solidearth.surfaceload.icethicknesschange(posant,:)=-50;

%greenland
posgre=find(late > 70 &  late < 80 & longe>-50 & longe<-30);
%md.solidearth.surfaceload.icethicknesschange(posgre,:)=-100;

%alaska : 
posala=find(late > 55 &  late < 68 & longe>-162 & longe<-140);
%md.solidearth.surfaceload.icethicknesschange(posala,:)=-150;


%elastic loading from love numbers:
md.solidearth.lovenumbers=lovenumbers('maxdeg',100);

%mask
mask=gmtmask(md.mesh.lat,md.mesh.long);
icemask=ones(md.mesh.numberofvertices,1);
pos=find(mask==0);  icemask(pos)=-1;
pos=find(sum(mask(md.mesh.elements),2)<3);   icemask(md.mesh.elements(pos,:))=-1;
md.mask.ice_levelset=icemask;
md.mask.ocean_levelset=-icemask;

%make sure that the elements that have loads are fully grounded:
%pos=find(sum(md.solidearth.surfaceload.icethicknesschange(1:end-1,:),2));
%md.mask.ocean_levelset(md.mesh.elements(pos,:))=1;

%make sure wherever there is an ice load, that the mask is set to ice:
md.mask.ice_levelset(md.mesh.elements(pos,:))=-1;

%dsl
md.dsl.global_average_thermosteric_sea_level_change=[0;0];
md.dsl.sea_surface_height_change_above_geoid=zeros(md.mesh.numberofvertices+1,1);
md.dsl.sea_water_pressure_change_at_sea_floor=zeros(md.mesh.numberofvertices+1,1);

md.solidearth.settings.ocean_area_scaling=0;

%Geometry for the bed, arbitrary: 
md.geometry.bed=-ones(md.mesh.numberofvertices,1);

%Materials: 
md.materials=materials('hydro');

%Miscellaneous
md.miscellaneous.name='scale2002';

%Uncertainty Quantification%
md.qmu.variables=struct();;

%partitioning scaling
npart=3;
partition=-ones(md.mesh.numberofelements,1);
partition(posgre)=0;
partition(posant)=1;
partition(posala)=2;

%variable scaling
md.qmu.variables.surfaceload0=normal_uncertain(...
	'descriptor','scaled_SurfaceloadIceThicknessChange',...
	'mean',ones(npart,nsteps),...
	'stddev',.1*ones(npart,nsteps),...
	'partition',partition,'nsteps',nsteps);

%responses 
md.qmu.responses.sealevel1=response_function('descriptor','Outputdefinition1');

%output definitions: 
md.outputdefinition.definitions={nodalvalue('name','SNode','definitionstring','Outputdefinition1', ...
			'model_string','Sealevel','node',1)};

%algorithm: 
md.qmu.method     =dakota_method('nond_samp');
md.qmu.method(end)=dmeth_params_set(md.qmu.method(end),...
	'seed',1234,...
	'samples',3,...
	'sample_type','random');

%parameters
md.qmu.params.direct=true;
md.qmu.params.interval_type='forward';
md.qmu.params.analysis_driver='matlab';
md.qmu.params.evaluation_scheduling='master';
md.qmu.params.processors_per_evaluation=1;
md.qmu.params.tabular_graphics_data=true;
md.qmu.output=1;

%transient: 
md.solidearth.settings.computesealevelchange=1;
md.timestepping.start_time=1999.5;
md.timestepping.interp_forcings=0;
md.timestepping.time_step=1;
md.timestepping.final_time=2099.5
md.transient.issmb=0;
md.transient.ismasstransport=0;
md.transient.isstressbalance=0;
md.transient.isthermal=0;
md.transient.isslr=1;
md.transient.isgia=1;
md.solidearth.requested_outputs= {'default',...
		'SurfaceloadIceThicknessChange','Sealevel','SealevelRSLRate',...
		'SealevelNEsaRate', 'SealevelUEsaRate', 'NGiaRate', 'UGiaRate',...
		'SealevelEustaticMask','SealevelEustaticOceanMask','MaskOceanLevelset','MaskIceLevelset'};

%hack: 
md.geometry.thickness=ones(md.mesh.numberofvertices,1);
md.geometry.base=-ones(md.mesh.numberofvertices,1);
md.geometry.surface=zeros(md.mesh.numberofvertices,1);

%settings: 
md.verbose.qmu=1;
md.verbose=verbose(0);
md.qmu.isdakota=1;

md.solidearth.settings.rigid=1; md.solidearth.settings.elastic=1;md.solidearth.settings.rotation=1;

md=solve(md,'tr');

%check scaling worked ok: 
n=md.qmu.method.params.samples;
field_values={};

for i=1:n,
	md2=model(); md2.results=md.results.dakota.modelresults{i};
	uq=md2.results.TransientSolution(1).uq_variables;
	part1=uq(1:100);
	part2=uq(101:200);
	part3=uq(201:300);

	h2=ones(md.mesh.numberofelements+1,100);
	
	for j=1:length(posgre),
		for k=1:100,
			h2(posgre(j),k)= h2(posgre(j),k)*part1(k);
		end
	end
	for j=1:length(posant),
		for k=1:100,
			h2(posant(j),k)= h2(posant(j),k)*part2(k);
		end
	end
	for j=1:length(posala),
		for k=1:100,
			h2(posala(j),k)= h2(posala(j),k)*part3(k);
		end
	end

	h=resultstomatrix(md2,'TransientSolution','SurfaceloadIceThicknessChange');
	h2(end,:)=h(end,:);

	field_values{end+1}=size(find(h2-h));
end


%Fields and tolerances to track changes
field_names     ={'dH1','dH2','dH3'};
field_tolerances={1e-5,1e-5,1e-5};

