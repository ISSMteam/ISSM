%Test Name: EarthSlr

%mesh earth:
md=model;
md.mesh=gmshplanet('radius',6.371012*10^3,'resolution',1500.); %500 km resolution mesh
md.cluster.np=16;

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
md.slr.hydro_rate = zeros(md.mesh.numberofvertices,1);

%GIA: 
md.gia=giamme();
md.gia.Ngia= ones(md.mesh.numberofvertices,10);
md.gia.Ugia= ones(md.mesh.numberofvertices,10);
load post;
ns=length(p);
for i=1:ns,
	md.gia.Ngia(:,i)=i;
	md.gia.Ugia(:,i)=i;
end
md.gia.modelid=10;
md.verbose=verbose(0);
md.verbose.qmu=1;

md.slr.requested_outputs= {'default',...
		'SealevelriseDeltathickness','Sealevel','SealevelRSLRate','SealevelriseCumDeltathickness',...
		'SealevelNEsaRate', 'SealevelUEsaRate', 'NGiaRate', 'UGiaRate',...
		'SealevelEustaticMask','SealevelEustaticOceanMask'};

ids=(1:(ns+1))';
probs=[p; 0]; 
%probs=[(1:ns)'; 0];  probs=probs/sum(probs);
md.qmu.variables.giamodelid=histogram_bin_uncertain('GiaModelid',ns+1,ids,probs);

%qmu: 
%x=0:.1:100;
%y=lognormal_pdf(x2,0,1);
%md.qmu.variables.surface_mass_balance=histogram_bin_uncertain('scaled_SmbMassBalance',length(x),x,[y 0]);

%x=1:(ns+1);
%y=[lognormal_pdf(x(1:end-1)+.5,0,1) 0]; y=y/sum(y);
%md.qmu.variables.giamodelid=histogram_bin_uncertain('GiaModelid',length(x),x,y);

%md.qmu.variables.giamodelid=uniform_uncertain('descriptor','GiaModelid','lower',1,...
%																		'upper',10);

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
'samples',50,...
'sample_type','lhs');

%parameters
md.qmu.params.direct=true;
md.qmu.params.interval_type='forward';
md.qmu.params.analysis_driver='matlab';
md.qmu.params.evaluation_scheduling='master';
md.qmu.params.processors_per_evaluation=2;
md.qmu.params.tabular_graphics_data=true;

md.qmu.isdakota=1;
md.qmu.output=1;

%Miscellaneous
md.miscellaneous.name='test2002';

%Solution parameters
md.slr.reltol=NaN;
md.slr.abstol=1e-3;
md.slr.geodetic=1;

% max number of iteration reverted back to 10 (i.e., the original default value)
md.slr.maxiter=10;

%eustatic run:
md.slr.rigid=1; md.slr.elastic=1;md.slr.rotation=1;
md=solve(md,'Sealevelrise');

