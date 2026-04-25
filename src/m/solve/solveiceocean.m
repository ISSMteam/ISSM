function md=solveiceocean(md,solutionstring,varargin)
%SOLVE - apply ice/ocean solution sequence for this model
%
%   Usage:
%      md=solveiceocean(md,solutionstring,varargin)
%      where varargin is a lit of paired arguments of string OR enums
%
%   solution types available comprise:
%		 - 'Stressbalance'      or 'sb'
%		 - 'Masstransport'      or 'mt'
%		 - 'Transient'          or 'tr'
%
%  extra options:
%      - loadonly    : does not solve. only load results
%      - runtimename : true or false (default is true), makes name unique
%      - checkconsistency : 'yes' or 'no' (default is 'yes'), ensures checks on consistency of model
%      - restart: 'directory name (relative to the execution directory) where the restart file is located.
%      - outbinread  : if 0, download the outbin but do not process is (md.results is not updated)
%
%   Examples:
%      md=solve(md,'Transient');
%      md=solve(md,'tr');

if ~ischar(solutionstring)
	error('ISSM''s solve function only accepts strings for solution sequences. Type help solve to get a list of supported solutions');
end

%recover and process solve options
if strcmpi(solutionstring,'sb') || strcmpi(solutionstring,'Stressbalance')
	solutionstring = 'StressbalanceSolution';
elseif strcmpi(solutionstring,'mt') || strcmpi(solutionstring,'Masstransport')
	solutionstring = 'MasstransportSolution';
elseif strcmpi(solutionstring,'tr') || strcmpi(solutionstring,'Transient')
	solutionstring = 'TransientSolution';
else
	error(['solutionstring ' solutionstring ' not supported!']);
end
options=pairoptions(varargin{:},'solutionstring',solutionstring);

%recover some fields
md.private.solution=solutionstring;
cluster=md.cluster;
if strcmpi(getfieldvalue(options,'batch','no'),'yes') batch=1; else batch=0; end
if ~isa(cluster,'localpfe') & ~isa(cluster,'generic')
	error('cannot run ice/ocean simulation on any other cluster than localpfe');
end

%check model consistency
if strcmpi(getfieldvalue(options,'checkconsistency','yes'),'yes')
	if md.verbose.solution,
		disp('checking model consistency');
	end
	ismodelselfconsistent(md),
end

%If we are restarting, actually use the provided runtime name:
restart=getfieldvalue(options,'restart','');
%First, build a runtime name that is unique
if restart==1 
	%Leave the runtimename as is
else
	if ~isempty(restart)
		md.private.runtimename=restart;
	elseif getfieldvalue(options,'runtimename',true)
		c=clock;
		md.private.runtimename=sprintf('%s-%02i-%02i-%04i-%02i-%02i-%02i-%i',md.miscellaneous.name,c(2),c(3),c(1),c(4),c(5),floor(c(6)),feature('GetPid'));
	else
		md.private.runtimename=md.miscellaneous.name;
	end
end

%if running qmu analysis, some preprocessing of dakota files using models
%fields needs to be carried out. 
if md.qmu.isdakota
	md=preqmu(md,options);
end

%Do we load results only?
if getfieldvalue(options,'loadonly',false),
	md=loadresultsfromcluster(md);
	return;
end

%Write all input files
marshall(md);                                          % bin file
ToolkitsFile(md.toolkits,[md.miscellaneous.name '.toolkits']); % toolkits file
BuildQueueScriptIceOcean(cluster,md.private.runtimename,md.miscellaneous.name,md.private.solution,md.settings.io_gather,md.debug.valgrind,md.debug.gprof,md.qmu.isdakota); % queue file

%Upload all required files
modelname = md.miscellaneous.name;
filelist  = {[modelname '.bin'] [modelname '.toolkits']};
if ispc
	filelist{end+1}=[modelname '.bat'];
else
	filelist{end+1}=[modelname '.queue'];
end

if md.qmu.isdakota,
	filelist{end+1} = [modelname '.qmu.in'];
end

if isempty(restart)
	UploadQueueJob(cluster,md.miscellaneous.name,md.private.runtimename,filelist);
end

%launch queue job: 
disp('launching solution sequence')
LaunchQueueJobIceOcean(cluster,md.miscellaneous.name,md.private.runtimename,filelist,restart,batch);

%return if batch: 
if batch,
	if md.verbose.solution,
		disp('batch mode requested: not launching job interactively');
		disp('launch solution sequence on remote cluster by hand');
	end
	return;
end
%wait on lock
if isnan(md.settings.waitonlock)
	%load when user enters 'y'
	disp('solution launched on remote cluster. log in to detect job completion.');
	choice=input('Is the job successfully completed? (y/n)','s');
	if ~strcmp(choice,'y'), 
		disp('Results not loaded... exiting'); 
	else
		md=loadresultsfromcluster(md);
	end
elseif md.settings.waitonlock>0
	%we wait for the done file
	done=waitonlock(md);
	if md.verbose.solution,
		disp('loading results from cluster');
	end
	md=loadresultsfromcluster(md,'runtimename','');
elseif md.settings.waitonlock==0
	 disp('Model results must be loaded manually with md=loadresultsfromcluster(md);');
end
