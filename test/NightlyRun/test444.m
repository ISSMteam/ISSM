%Test Name: SquareSheetShelfTranSSA2dAggressiveDakotaSampRegionalOutput

% TODO:
% - Figure out why test fails intermittently on Mac with "Index exceeds array 
% bounds."
%

md=triangle(model(),'../Exp/Square.exp',150000.);
md=setmask(md,'../Exp/SquareShelf.exp','');
md=parameterize(md,'../Par/SquareSheetShelf.par');
md.geometry.bed=md.geometry.base;
pos=find(md.mask.ocean_levelset<0);
md.geometry.bed(pos)=md.geometry.base(pos)-10;
md.friction.coefficient=20.*ones(md.mesh.numberofvertices,1);
md.friction.p=ones(md.mesh.numberofelements,1);
md.friction.q=ones(md.mesh.numberofelements,1);
md.transient.isthermal=0;
md.transient.isgroundingline=1;
md.groundingline.migration='AggressiveMigration';

md.settings.output_frequency=3;
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',3);

regionalmask=zeros(md.mesh.numberofvertices,1);
in=ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,'../Exp/SquareHalfRight.exp','node',1);
regionalmask(find(in))=1;
md.transient.requested_outputs={'default','GroundedArea','FloatingArea','IceVolumeAboveFloatation','GroundedArea1','FloatingArea1','TotalFloatingBmb1','TotalGroundedBmb1','TotalSmb1',...
	'IceMass1','IceVolume1','IceVolumeAboveFloatation1','IceVolumeAboveFloatation'};
md.outputdefinition.definitions{1}=regionaloutput('name','GroundedArea1','outputnamestring','GroundedArea','mask',regionalmask,...
	'definitionstring','Outputdefinition1');
md.outputdefinition.definitions{2}=regionaloutput('name','FloatingArea1','outputnamestring','FloatingArea','mask',regionalmask,...
	'definitionstring','Outputdefinition2');
md.outputdefinition.definitions{3}=regionaloutput('name','TotalFloatingBmb1','outputnamestring','TotalFloatingBmb','mask',regionalmask,...
	'definitionstring','Outputdefinition3');
md.outputdefinition.definitions{4}=regionaloutput('name','TotalGroundedBmb1','outputnamestring','TotalGroundedBmb','mask',regionalmask,...
	'definitionstring','Outputdefinition4');
md.outputdefinition.definitions{5}=regionaloutput('name','IceMass1','outputnamestring','IceMass','mask',regionalmask,...
	'definitionstring','Outputdefinition5');
md.outputdefinition.definitions{6}=regionaloutput('name','IceVolume1','outputnamestring','IceVolume','mask',regionalmask,...
	'definitionstring','Outputdefinition6');
md.outputdefinition.definitions{7}=regionaloutput('name','IceVolumeAboveFloatation1','outputnamestring','IceVolumeAboveFloatation','mask',regionalmask,...
	'definitionstring','Outputdefinition7');
md.outputdefinition.definitions{8}=regionaloutput('name','TotalSmb1','outputnamestring','TotalSmb','mask',regionalmask,...
	'definitionstring','Outputdefinition8');
md.outputdefinition.definitions{9}=regionaloutput('name','TotalSmb2','outputnamestring','TotalSmb','mask',regionalmask,...
	 'definitionstring','Outputdefinition9');

md=extrude(md,3,1);
md=collapse(md);

%Dakota options

%dakota version
version=IssmConfig('_DAKOTA_VERSION_'); version=version(1:3); version=str2num(version);

%partitioning
npart=10;
partition=partitioner(md,'package','chaco','npart',npart,'weighting','on')-1;
md.qmu.isdakota=1;

md.qmu.variables.drag_coefficient=normal_uncertain('descriptor','scaled_BasalforcingsFloatingiceMeltingRate',...
	'mean',ones(npart,1),...
	'stddev',.1*ones(npart,1),...
	'partition',partition);

md.qmu.responses.IceMass1=response_function('descriptor','Outputdefinition5');
md.qmu.responses.IceVolume1=response_function('descriptor','Outputdefinition6');
md.qmu.responses.IceVolumeAboveFloatation1=response_function('descriptor','Outputdefinition7');
md.qmu.responses.IceVolumeAboveFloatation=response_function('descriptor','IceVolumeAboveFloatation');
md.qmu.responses.GroundedArea1=response_function('descriptor','Outputdefinition1');
md.qmu.responses.FloatingArea1=response_function('descriptor','Outputdefinition2');
md.qmu.responses.TotalFloatingBmb1=response_function('descriptor','Outputdefinition3');
md.qmu.responses.TotalGroundedBmb1=response_function('descriptor','Outputdefinition4');
md.qmu.responses.TotalSmb1=response_function('descriptor','Outputdefinition8');
md.qmu.responses.TotalSmb2=response_function('descriptor','Outputdefinition9');
md.qmu.responses.FloatingArea=response_function('descriptor','FloatingArea');

md.qmu.method     =dakota_method('nond_samp');
md.qmu.method(end)=dmeth_params_set(md.qmu.method(end),...
	'seed',1234,...
	'samples',20,...
	'sample_type','random');

%%  a variety of parameters
md.qmu.params.direct=true;
md.qmu.params.analysis_components='';
md.qmu.params.tabular_graphics_data=true;

if version>=6,
	md.qmu.params.analysis_driver='matlab';
	md.qmu.params.evaluation_scheduling='master';
	md.qmu.params.processors_per_evaluation=2;
else
	md.qmu.params.analysis_driver='stressbalance';
	md.qmu.params.evaluation_concurrency=1;
end


md.stressbalance.reltol=10^-5; %tighten for qmu analyses

md=solve(md,'Transient','overwrite','y');

%Fields and tolerances to track changes
md.qmu.results=md.results.dakota;

%we put all the mean and stdev data in the montecarlo field, which we will use to test for success.
md.results.dakota.montecarlo=[];
for i=1:11,
	md.results.dakota.montecarlo=[md.results.dakota.montecarlo md.results.dakota.dresp_out(i).mean];
end
for i=1:11,
	md.results.dakota.montecarlo=[md.results.dakota.montecarlo md.results.dakota.dresp_out(i).stddev];
end
field_names     ={'montecarlo'};
field_tolerances={1e-11};
field_values={...
	md.results.dakota.montecarlo,...
	};

