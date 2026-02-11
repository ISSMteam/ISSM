function slm=solveslm(slm,solutionstringi,varargin)
%SOLVESLM - apply solution sequence for this sealevel model
%
%   Usage:
%      slm=solveslm(slm,solutionstring,varargin)
%      where varargin is a lit of paired arguments of string OR enums
%
%   solution types available comprise:
%      - 'Transient'
%
%   extra options:
%
%   Examples:
%      slm=solveslm(slm,'Transient');

%recover and process solve options
if strcmpi(solutionstringi,'tr') || strcmpi(solutionstringi,'Transient')
	solutionstring = 'TransientSolution';
else
	error(['solutionstring ' solutionstringi ' not supported!']);
end

%Default settings for debugging
valgrind=0;
%slm.cluster.interactive=0; valgrind=1;

%check consistency: 
slm.checkconsistency(solutionstring);

%process options 
options=pairoptions(varargin{:},'solutionstring',solutionstring);

%make sure we request sum of cluster processors 
totalnp=0;
for i=1:length(slm.icecaps), totalnp=totalnp+slm.icecaps{i}.cluster.np; end
totalnp=totalnp+slm.earth.cluster.np;
if totalnp~=slm.cluster.np,
	error('sum of all icecaps and earch cluster processors requestes should be equal to slm.cluster.np');
end

%recover some fields
slm.private.solution=solutionstring;
cluster=slm.cluster;
batch=0;
%now, go through icecaps, glaciers and earth, and upload all the data independently: 
disp('solving ice caps first');
for i=1:length(slm.icecaps),
	slm.icecaps{i}=solve(slm.icecaps{i},solutionstringi,'batch','yes');
end
disp('solving earth now');
slm.earth=solve(slm.earth,solutionstringi,'batch','yes');

%First, build a runtime name that is unique
c=clock;
slm.private.runtimename=sprintf('%s-%02i-%02i-%04i-%02i-%02i-%02i-%i',slm.miscellaneous.name,c(2),c(3),c(1),c(4),c(5),floor(c(6)),feature('GetPid'));

%Write all input files:
privateruntimenames={}; 
miscellaneousnames={}; 
nps={};
for i=1:length(slm.icecaps),
	privateruntimenames{end+1}=slm.icecaps{i}.private.runtimename;
	miscellaneousnames{end+1}=slm.icecaps{i}.miscellaneous.name;
	nps{end+1}=slm.icecaps{i}.cluster.np;
end
privateruntimenames{end+1}=slm.earth.private.runtimename;
miscellaneousnames{end+1}=slm.earth.miscellaneous.name;
nps{end+1}=slm.earth.cluster.np;

BuildQueueScriptMultipleModels(cluster,slm.private.runtimename,slm.miscellaneous.name,slm.private.solution,privateruntimenames,miscellaneousnames,nps);

%Upload all required files, given that each individual solution for icecaps and earth model already did:
filelist={[slm.miscellaneous.name '.queue']};
UploadQueueJob(cluster,slm.miscellaneous.name,slm.private.runtimename,filelist);

%launch queue job: 
disp('launching solution sequence')
LaunchQueueJob(cluster,slm.miscellaneous.name,slm.private.runtimename,filelist,'',batch);

%wait on lock
if isnan(slm.settings.waitonlock),
	%load when user enters 'y'
	disp('solution launched on remote cluster. log in to detect job completion.');
	choice=input('Is the job successfully completed? (y/n)','s');
	if ~strcmp(choice,'y'), 
		disp('Results not loaded... exiting'); 
	else
		for i=1:length(slm.icecaps), slm.icecaps{i}=loadresultsfromcluster(slm.icecaps{i});end;
		slm.earth=loadresultsfromcluster(slm.earth);
	end
elseif slm.settings.waitonlock>0,
	%we wait for the done file
	done=waitonlock(slm);
	disp('loading results from cluster');
	for i=1:length(slm.icecaps), slm.icecaps{i}=loadresultsfromcluster(slm.icecaps{i});end;
	slm.earth=loadresultsfromcluster(slm.earth);
elseif slm.settings.waitonlock==0,
	 disp('Model results must be loaded manually with slm=loadresultsfromcluster(slm);');
end
