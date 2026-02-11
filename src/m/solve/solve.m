function md=solve(md,solutionstring,varargin)
%SOLVE - apply solution sequence for this model
%
%   Usage:
%      md=solve(md,solutionstring,varargin)
%
%   where varargin is a lit of paired arguments of string OR enums
%
%   Solution types available comprise:
%   - 'Stressbalance'        or 'sb'
%   - 'Masstransport'        or 'mt'
%   - 'Oceantransport'       or 'oceant'
%   - 'Thermal'              or 'th'
%   - 'Steadystate'          or 'ss'
%   - 'Transient'            or 'tr'
%   - 'Balancethickness'     or 'mc'
%   - 'Balancethickness2'
%   - 'BalancethicknessSoft' or 'mcsoft'
%   - 'Balancevelocity'      or 'bv'
%   - 'BedSlope'             or 'bsl'
%   - 'SurfaceSlope'         or 'ssl'
%   - 'Hydrology'            or 'hy'
%   - 'DamageEvolution'      or 'da'
%   - 'Gia'                  or 'gia'
%   - 'Love'                 or 'lv'
%   - 'Esa'                  or 'esa'
%   - 'Sampling'             or 'smp'
%
%   Extra options:
%   - loadonly         : do not solve, only load results
%   - runtimename      : true or false (default is true); makes name unique
%   - checkconsistency : 'yes' or 'no' (default is 'yes'); checks consistency of model
%   - restart          : directory name (relative to the execution directory) 
%                        where the restart file is located
%
%   Examples:
%      md=solve(md,'Stressbalance');
%      md=solve(md,'sb');

if ~ischar(solutionstring)
	error('ISSM''s solve function only accepts strings for solution sequences. Type help solve to get a list of supported solutions');
end

%recover and process solve options
if strcmpi(solutionstring,'sb') || strcmpi(solutionstring,'Stressbalance')
	solutionstring = 'StressbalanceSolution';
elseif strcmpi(solutionstring,'mt') || strcmpi(solutionstring,'Masstransport')
	solutionstring = 'MasstransportSolution';
elseif strcmpi(solutionstring,'oceant') || strcmpi(solutionstring,'Oceantransport')
	solutionstring = 'OceantransportSolution';
elseif strcmpi(solutionstring,'th') || strcmpi(solutionstring,'Thermal')
	solutionstring = 'ThermalSolution';
elseif strcmpi(solutionstring,'ss') || strcmpi(solutionstring,'Steadystate') 
	solutionstring = 'SteadystateSolution';
elseif strcmpi(solutionstring,'tr') || strcmpi(solutionstring,'Transient')
	solutionstring = 'TransientSolution';
elseif strcmpi(solutionstring,'mc') || strcmpi(solutionstring,'Balancethickness')
	solutionstring = 'BalancethicknessSolution';
elseif strcmpi(solutionstring,'Balancethickness2')
	solutionstring = 'Balancethickness2Solution';
elseif strcmpi(solutionstring,'mcsoft') || strcmpi(solutionstring,'BalancethicknessSoft')
	solutionstring = 'BalancethicknessSoftSolution';
elseif strcmpi(solutionstring,'bv') || strcmpi(solutionstring,'Balancevelocity')
	solutionstring = 'BalancevelocitySolution';
elseif strcmpi(solutionstring,'bsl') || strcmpi(solutionstring,'BedSlope')
	solutionstring = 'BedSlopeSolution';
elseif strcmpi(solutionstring,'ssl') || strcmpi(solutionstring,'SurfaceSlope')
	solutionstring = 'SurfaceSlopeSolution';
elseif strcmpi(solutionstring,'hy') || strcmpi(solutionstring,'Hydrology')
	solutionstring = 'HydrologySolution';
elseif strcmpi(solutionstring,'da') || strcmpi(solutionstring,'DamageEvolution')
	solutionstring = 'DamageEvolutionSolution';
elseif strcmpi(solutionstring,'gia') || strcmpi(solutionstring,'Gia')
	solutionstring = 'GiaSolution';
elseif strcmpi(solutionstring,'lv') || strcmpi(solutionstring,'Love')
	solutionstring = 'LoveSolution';
elseif strcmpi(solutionstring,'esa') || strcmpi(solutionstring,'Esa')
	solutionstring = 'EsaSolution';
elseif strcmpi(solutionstring,'smp') || strcmpi(solutionstring,'Sampling')
	solutionstring = 'SamplingSolution';    
else
	error(['solutionstring ' solutionstring ' not supported!']);
end
options=pairoptions(varargin{:},'solutionstring',solutionstring);

%recover some fields
md.private.solution=solutionstring;
cluster=md.cluster;
if strcmpi(getfieldvalue(options,'batch','no'),'yes')
	batch=1;
else
	batch=0;
end

%check model consistency
if strcmpi(getfieldvalue(options,'checkconsistency','yes'),'yes'),
	if md.verbose.solution,
		disp('checking model consistency');
	end
	ismodelselfconsistent(md);
end

%If we are restarting, actually use the provided runtime name:
restart=getfieldvalue(options,'restart','');
%First, build a runtime name that is unique
if restart==1 
	%Leave the runtimename as is
else
	if ~isempty(restart),
		md.private.runtimename=restart;
	elseif getfieldvalue(options,'runtimename',true),
		c=clock;
		md.private.runtimename=sprintf('%s-%02i-%02i-%04i-%02i-%02i-%02i-%i',md.miscellaneous.name,c(2),c(3),c(1),c(4),c(5),floor(c(6)),feature('GetPid'));
	else
		md.private.runtimename=md.miscellaneous.name;
	end
end

%if running QMU analysis, some preprocessing of Dakota files using model fields needs to be carried out. 
if md.qmu.isdakota,
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
BuildQueueScript(cluster,md.private.runtimename,md.miscellaneous.name,md.private.solution,md.settings.io_gather,md.debug.valgrind,md.debug.gprof,md.qmu.isdakota,md.transient.isoceancoupling); % queue file

%Upload all required files
modelname = md.miscellaneous.name;
filelist  = {[modelname '.bin'] [modelname '.toolkits']};
if ispc,
	filelist{end+1}=[modelname '.bat'];
else
	filelist{end+1}=[modelname '.queue'];
end

if md.qmu.isdakota,
	filelist{end+1} = [modelname '.qmu.in'];
end

if isempty(restart),
	disp('uploading input files')
	UploadQueueJob(cluster,md.miscellaneous.name,md.private.runtimename,filelist);
end

%launch queue job: 
disp('launching solution sequence')
LaunchQueueJob(cluster,md.miscellaneous.name,md.private.runtimename,filelist,restart,batch);

%return if batch: 
if batch
	if md.verbose.solution
		disp('batch mode requested: not launching job interactively');
		disp('launch solution sequence on remote cluster by hand');
	end
	return;
end

%wait on lock
if isnan(md.settings.waitonlock),
	%load when user enters 'y'
	disp('solution launched on remote cluster. log in to detect job completion.');
	choice=input('Is the job successfully completed? (y/n)','s');
	if ~strcmp(choice,'y'), 
		disp('Results not loaded... exiting'); 
	else
		md=loadresultsfromcluster(md);
	end
elseif md.settings.waitonlock>0,
	%wait for done file
	done=waitonlock(md);
	if md.verbose.solution,
		disp('loading results from cluster');
	end
	md=loadresultsfromcluster(md);
elseif md.settings.waitonlock==0,
	 disp('Model results must be loaded manually with md=loadresultsfromcluster(md);');
end
