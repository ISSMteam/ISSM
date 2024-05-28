function ispresent=waitonlock(md)
%WAITONLOCK - wait for a file
%
%   This routine will return when a file named 'lockfilename' is written to 
%   disk. Also check for outlog file be cause it might be written several 
%   seconds after the lock file.
%
%   If the time limit given in input is exceeded, return 0
%
%   Usage:
%      flag=waitonlock(md)

%Return if waitonlock < 0 (no need to wait)

%Get lockfilename (lock file) and options
executionpath = md.cluster.executionpath;
timelimit     = md.settings.waitonlock;
cluster       = md.cluster;

if isa(cluster,'pfe') && cluster.interactive>0
	lockfilename  = [executionpath '/Interactive' num2str(cluster.interactive) '/' md.miscellaneous.name '.lock'];
	logfilename   = [executionpath '/Interactive' num2str(cluster.interactive) '/' md.miscellaneous.name '.outlog'];
elseif isa(cluster,'localpfe'),
	lockfilename  = [executionpath '/' md.miscellaneous.name '.lock'];
	logfilename   = [executionpath '/' md.miscellaneous.name '.outlog'];
else
	lockfilename  = [executionpath '/' md.private.runtimename '/' md.miscellaneous.name '.lock'];
	logfilename   = [executionpath '/' md.private.runtimename '/' md.miscellaneous.name '.outlog'];
end


%If we are using the generic cluster in interactive mode, job is already complete
if (isa(cluster,'generic') & cluster.interactive) | isa(cluster,'generic_static'),
	%We are in interactive mode, no need to check for job completion
	ispresent=1;
	return;
end

%initialize time and file presence test flag
elapsedtime=0; ispresent=0; starttime=clock;
disp(['waiting for ' lockfilename ' hold on... (Ctrl+C to exit)'])

%prepare command if the job is not running on the local machine
if ~strcmpi(oshostname(),cluster.name),
	if isa(cluster,'cloud'),
		command = [' [ -f ' lockfilename ' ] && [ -f ' logfilename ' ] 2>/dev/null'];
		command = [starcluster() ' sshmaster ' cluster.name ' --user ' cluster.login ' ''' command ''''];
	else
		command = ['ssh -l ' cluster.login];
		if isprop(cluster,'idfile') && ~strcmp(cluster.idfile,''),
			command = [command ' -i ' cluster.idfile];
		end
		if isprop(cluster,'port') && cluster.port,
			command = [command ' -p ' num2str(cluster.port) ' localhost'];
		else,
			command = [command ' ' cluster.name];
		end
		command = [command ' "[ -f ' lockfilename ' ] && [ -f ' logfilename ' ]" 2>/dev/null'];
	end
end

%loop till file .lock exist or time is up
while (ispresent==0 & elapsedtime<timelimit)
	if strcmpi(oshostname(),cluster.name),
		pause(1);
		ispresent=(exist(lockfilename,'file') & exist(logfilename,'file'));
		elapsedtime=etime(clock,starttime)/60;
	else
		pause(5);
		elapsedtime=etime(clock,starttime);
		fprintf('\rchecking for job completion (time: %i min %i sec)      ',floor(elapsedtime/60),floor(rem(elapsedtime,60)));
		elapsedtime=elapsedtime/60; %converts time from sec to min
		ispresent=~system(command);
		if ispresent, fprintf('\n'); end
	end
end

%build output
if (elapsedtime>timelimit),
	disp('Time limit exceeded. Increase md.settings.waitonlock');
	disp('The results must be loaded manually with md=loadresultsfromcluster(md).');
	error(['waitonlock error message: time limit exceeded']);
end
