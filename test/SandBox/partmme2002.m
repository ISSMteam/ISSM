%Test Name: EarthSlr

%mesh earth:
md=model;
md.mesh=gmshplanet('radius',6.371012*10^3,'resolution',500.); %500 km resolution mesh
md.cluster.np=2;

%parameterize slr solution:
%slr loading:  {{{
nsteps=100;
md.solidearth.surfaceload.icethicknesschange=zeros(md.mesh.numberofelements+1,nsteps);
md.solidearth.surfaceload.icethicknesschange(end,:)=2000.5:1:2099.5;
md.solidearth.sealevel=zeros(md.mesh.numberofvertices,1);
md.dsl.global_average_thermosteric_sea_level_change=[0;0];
md.dsl.sea_surface_height_change_above_geoid=zeros(md.mesh.numberofvertices+1,1);
md.dsl.sea_water_pressure_change_at_sea_floor=zeros(md.mesh.numberofvertices+1,1);

%antarctica
late=sum(md.mesh.lat(md.mesh.elements),2)/3;
longe=sum(md.mesh.long(md.mesh.elements),2)/3;

posant=find(late <-80);
md.solidearth.surfaceload.icethicknesschange(posant,:)=-50;

%greenland
posgre=find(late > 70 &  late < 80 & longe>-60 & longe<-30);
md.solidearth.surfaceload.icethicknesschange(posgre,:)=-100;

%alaska : 
posala=find(late > 62 &  late < 68 & longe>-162 & longe<-140);
md.solidearth.surfaceload.icethicknesschange(posala,:)=-150;


%hawaii : 
poshaw=find(late > 15 &  late < 25 & longe>-170 & longe<-130);
md.solidearth.surfaceload.icethicknesschange(poshaw,:)=10;

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
pos=find(sum(md.solidearth.surfaceload.icethicknesschange(1:end-1,:),2));
md.mask.ocean_levelset(md.mesh.elements(pos,:))=1;

%make sure wherever there is an ice load, that the mask is set to ice:
md.mask.ice_levelset(md.mesh.elements(pos,:))=-1;
% }}}

md.solidearth.settings.ocean_area_scaling=0;

%Geometry for the bed, arbitrary: 
md.geometry.bed=-ones(md.mesh.numberofvertices,1);

%Materials: 
md.materials=materials('hydro');

%Miscellaneous
md.miscellaneous.name='test2002';

%Solution parameters
%Solution parameters
md.solidearth.settings.reltol=NaN;
md.solidearth.settings.abstol=1e-3;
md.solidearth.settings.computesealevelchange=1;

% max number of iteration reverted back to 10 (i.e., the original default value)
md.solidearth.settings.maxiter=10;

%initialize GIA: 
md.gia=giamme();
md.gia.Ngia=zeros(md.mesh.numberofvertices,1);
md.gia.Ugia=zeros(md.mesh.numberofvertices,1);
md.gia.modelid=1;

%Uncertainty Quantification%{{{
md.qmu.variables=struct();;

%partition vector: 
npart=3;
partition=-ones(md.mesh.numberofelements,1);

partition(posant)=0;
partition(posgre)=1;
partition(posala)=2;

vposant=unique(md.mesh.elements(posant,:));
vposgre=unique(md.mesh.elements(posgre,:));
vposala=unique(md.mesh.elements(posala,:));


%prepare arrays: 
nant=10; ngre=20; nala=30; nmax=max([nant,ngre,nala]);

%make ocean levelset temporal:
md.mask.ocean_levelset=[md.mask.ocean_levelset;2000];

%background arrays: 
bkgd=md.solidearth.surfaceload.icethicknesschange; 
bkgdnan=bkgd; bkgdnan(1:end-1,:)=NaN;

bkgdocean=md.mask.ocean_levelset;
bkgdoceannan=bkgdocean; bkgdoceannan(1:end-1,1)=NaN;

%mme arrays:
md.solidearth.surfaceload.icethicknesschange=cell(30,1);
md.mask.ocean_levelset=cell(30,1);
for i=1:nmax,
	if i==1,
		md.solidearth.surfaceload.icethicknesschange{i}=bkgd;
		md.mask.ocean_levelset{i}=bkgdocean;
	else
		bi=bkgdnan;
		oi=bkgdoceannan;
		if i<=10, 
			bi(posant,:)=(-50+i)*ones(length(posant),nsteps); 
			%oi(vposant,:)=bkgdocean(vposant,:)-50+i;
			oi(vposant,:)=i;
		end
		if i<=20, 
			bi(posgre,:)=(-100+i)*ones(length(posgre),nsteps); 
			%oi(vposgre,:)=bkgdocean(vposgre,:)-100+i;
			oi(vposgre,:)=2*i;
		end
		if i<=30, 
			bi(posala,:)=(-150+i)*ones(length(posala),nsteps); 
			%oi(vposala,:)=bkgdocean(vposala,:)-150+i;
			oi(vposala,:)=3*i;
		end
		md.solidearth.surfaceload.icethicknesschange{i}=bi;
		ois=repmat(oi,1,nsteps); ois(end,:)=bi(end,:);
		md.mask.ocean_levelset{i}=ois;
	end
	md.mask.ocean_levelset{i}=bkgdocean;
end
%md.mask.ice_levelset=md.mask.ocean_levelset;
%md.solidearth.surfaceload.icethicknesschange=bkgd;
%md.mask.ocean_levelset=bkgdocean;

%create distributed histograms for each partition:
pairs_per_variable=zeros(npart,1);
abscissas=cell(npart,1);
counts=cell(npart,1);
[abscissas{1} counts{1} pairs_per_variable(1)]=equiprobable_histogram_uncertain(10); 
[abscissas{2} counts{2} pairs_per_variable(2)]=equiprobable_histogram_uncertain(20); 
[abscissas{3} counts{3} pairs_per_variable(3)]=equiprobable_histogram_uncertain(30); 

%variables: 
md.qmu.variables.surfaceload=histogram_bin_uncertain(...
	'descriptor','distributed_SurfaceloadModelid',...
	'pairs_per_variable',pairs_per_variable,'abscissas',abscissas,'counts',counts,'partition',partition);
md.qmu.correlation_matrix=[];

%responses 
md.qmu.responses.sealevel1=response_function('descriptor','Outputdefinition1');
md.qmu.responses.sealevel2=response_function('descriptor','Outputdefinition2');
md.qmu.responses.sealevel3=response_function('descriptor','Outputdefinition3');
md.qmu.responses.sealevel4=response_function('descriptor','Outputdefinition4');
md.qmu.responses.sealevel5=response_function('descriptor','Outputdefinition5');
md.qmu.responses.sealevel6=response_function('descriptor','Outputdefinition6');
md.qmu.responses.sealevel7=response_function('descriptor','Outputdefinition7');
md.qmu.responses.sealevel8=response_function('descriptor','Outputdefinition8');
md.qmu.responses.sealevel8=response_function('descriptor','Outputdefinition9');
md.qmu.responses.sealevel10=response_function('descriptor','Outputdefinition10');

%output definitions: 
for i=1:10,
	if i==1,
		md.outputdefinition.definitions={nodalvalue('name','SNode','definitionstring','Outputdefinition1', ...
			'model_string','Sealevel','node',i)}; 
	else
		md.outputdefinition.definitions{end+1}=nodalvalue('name','SNode','definitionstring',['Outputdefinition' num2str(i)], ...
			'model_string','Sealevel','node',i); 
end
end
%algorithm: 
md.qmu.method     =dakota_method('nond_samp');
md.qmu.method(end)=dmeth_params_set(md.qmu.method(end),...
	'seed',1234,...
	'samples',3,...
	'sample_type','lhs');

%parameters
md.qmu.params.direct=true;
md.qmu.params.interval_type='forward';
md.qmu.params.analysis_driver='matlab';
md.qmu.params.evaluation_scheduling='master';
md.qmu.params.processors_per_evaluation=1;
md.qmu.params.tabular_graphics_data=true;
md.qmu.output=1;
%}}}

%transient: 
md.timestepping.start_time=2000;
md.timestepping.interp_forcings=1;
md.timestepping.final_time=2002;
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

%eustatic + rigid + elastic + rotation run:
md.verbose=verbose('11111111111');
%md.verbose=verbose(0);
md.verbose.qmu=1;
md.solidearth.settings.rigid=1; md.solidearth.settings.elastic=1;md.solidearth.settings.rotation=1;
md.qmu.isdakota=1;
md=solve(md,'tr');

