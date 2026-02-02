%HPC class definition
%
%   Usage:
%      cluster=hpc_simba();
%      cluster=hpc_simba('np',3);
%      cluster=hpc_simba('np',3,'login','username');

classdef hpc_simba2
    properties (SetAccess=public)  
		 % {{{
		 name='simba20'
		 login='inwoo';
		 numnodes=18;    % number of nodes at 2019-11 installation
		 cpuspernode=36; % default number of cpus at each node
		 node=1;         % number of nodes for calculating 
		 np=4;           % number of cpus for calculating at each node
		 port=0;        % not know open port
		 queue='pub64';
		 codepath= [];
		 executionpath=[];
		 interactive=0;
		 verbose=0;     % show process of downloading
		 isqsub=1;
	 end
	 %}}}
	 methods
		 function cluster=hpc_simba2(varargin) % {{{

			 %initialize cluster using default settings if provided
			 if (exist('hpc_settings')==2), hpc_settings; end

			 %use provided options to change fields
			 cluster=AssignObjectFields(pairoptions(varargin{:}),cluster);

			 % where is "ISSM_DIR"?
			 if strcmpi(cluster.name,'simba00')
				 ISSM_DIR=getenv('ISSM_DIR');
			 elseif strcmpi(cluster.name,'simba20')
				 ISSM_DIR='/home/inwoo/issm/trunk-jpl/';
			 else
				 error(sprintf('ERROR: %s is not supported cluster name...'))
			 end
			 cluster.codepath=sprintf('%s/bin',ISSM_DIR);

			 % define specific user
			 [~, s] = system('whoami'); % get user name
			 s = s(1:end-1);
			 cluster.login=s;
			 if strcmpi(s,'inwoo'),
				 %cluster.executionpath=sprintf('%s/execution/',ISSM_DIR);
				 cluster.executionpath='/data2/inwoo/execution';
			 elseif strcmpi(s,'emilia'),
				 cluster.executionpath='/home/DATA/emilia-model/executionlog';
			 end
		 end
		 %}}}
		 function disp(cluster) % {{{
			 %  display the object
			 disp(sprintf('class ''%s'' object ''%s'' = ',class(cluster),inputname(1)));
			 disp(sprintf('    name: %s',cluster.name));
			 disp(sprintf('    login: %s',cluster.login));
			 disp(sprintf('    port: %i',cluster.port));
			 disp(sprintf('    numnodes: %i',cluster.numnodes));
			 disp(sprintf('    node: %i',cluster.node));
			 disp(sprintf('    cpuspernode: %i',cluster.cpuspernode));
			 disp(sprintf('    np: %i',cluster.np));
			 disp(sprintf('    queue: %s',cluster.queue));
			 disp(sprintf('    codepath: %s',cluster.codepath));
			 disp(sprintf('    executionpath: %s',cluster.executionpath));
			 disp(sprintf('    interactive: %i',cluster.interactive));
			 disp(sprintf('    isqub: %i',cluster.isqsub));
		 end
		 %}}}
		 function md = checkconsistency(cluster,md,solution,analyses) % {{{

			 available_queues={'pub64','free64','free48','free*,pub64','free*'};
			 queue_requirements_time=[Inf Inf Inf Inf Inf];
			 queue_requirements_np=[64 64 48 48 48];

			 QueueRequirements(available_queues,queue_requirements_time,queue_requirements_np,cluster.queue,cluster.np,1)

			 %Miscelaneous
			 if isempty(cluster.login), md = checkmessage(md,'login empty'); end
			 if isempty(cluster.codepath), md = checkmessage(md,'codepath empty'); end
			 if isempty(cluster.executionpath), md = checkmessage(md,'executionpath empty'); end

		 end
		 %}}}
		 function BuildKrigingQueueScript(cluster,modelname,solution,io_gather,isvalgrind,isgprof) % {{{

			 if(isvalgrind), disp('valgrind not supported by cluster, ignoring...'); end
			 if(isgprof),    disp('gprof not supported by cluster, ignoring...'); end

			 %write queuing script 
			 fid=fopen([modelname '.queue'],'w');
			 fprintf(fid,'#!/bin/bash\n');
			 fprintf(fid,'#$ -N %s\n',modelname);
			 fprintf(fid,'#$ -q %s \n',cluster.queue);
			 fprintf(fid,'#$ -pe one-node-mpi 2-64\n');
			 fprintf(fid,'#$ -R y\n');
			 fprintf(fid,'#$ -m beas\n');
			 fprintf(fid,'#$ -o %s.outlog \n',modelname);
			 fprintf(fid,'#$ -e %s.errlog \n\n',modelname);
			 fprintf(fid,'export ISSM_DIR="%s/../"\n',cluster.codepath); %FIXME
			 fprintf(fid,'source $ISSM_DIR/etc/environment.sh\n');       %FIXME
			 fprintf(fid,'cd %s/%s\n\n',cluster.executionpath,modelname);
			 fprintf(fid,'mpiexec -np %i %s/kriging.exe %s %s\n',cluster.np,cluster.codepath,[cluster.executionpath '/' modelname],modelname);
			 if ~io_gather, %concatenate the output files:
				 fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
			 end
			 fclose(fid);
		 end
		 %}}}
		 function BuildQueueScript(cluster,dirname,modelname,solution,io_gather,isvalgrind,isgprof,isdakota,isoceancoupling) % {{{

			 if(isvalgrind), disp('valgrind not supported by cluster, ignoring...'); end
			 if(isgprof),    disp('gprof not supported by cluster, ignoring...'); end

			 %write queuing script 
			 fid=fopen([modelname '.queue'],'w');
			 %fprintf(fid,'#PBS -q workq\n');
			 fprintf(fid,'#PBS -S /bin/bash\n');      % set shell name
			 fprintf(fid,'#PBS -N %s\n',modelname);   % set job name
			 % 8 node is available at simba
			 if (cluster.node >= cluster.numnodes),
				 fprintf('cluster node   : simba%02\n', cluster.node);
				 fprintf('number of node : %d\n', cluster.numnodes);
				 error('ERROR: cluster node is higher than number of nodes');
			 end
			 % node             -> how many nodes do you use?
			 % ppn (cluster.np) -> how many cpus do you use? 
			 fprintf(fid,'#PBS -l nodes=%d:ppn=%d\n',cluster.node,cluster.np);
			 fprintf(fid,'#PBS -o %s.outlog \n',modelname);
			 fprintf(fid,'#PBS -e %s.errlog \n',modelname);
			 fprintf(fid,'\n');
			 fprintf(fid,'source ~/.bashrc\n');                          % FIXME
			 fprintf(fid,'export ISSM_DIR="%s/../"\n',cluster.codepath); %FIXME
			 fprintf(fid,'source $ISSM_DIR/etc/environment.sh\n');       %FIXME
			 fprintf(fid,'\n');
			 fprintf(fid,'module load intel18/compiler-18\n');
			 fprintf(fid,'module load intel18/mvapich2-2.2\n');
			 fprintf(fid,'\n');
			 fprintf(fid,'cd %s/%s\n\n',cluster.executionpath,dirname);

			 if true, % HACK
				 % NOTE: old version
				 % fprintf(fid,'mpiexec -genv MV2_ENABLE_AFFINITY 0 -np %i %s/issm.exe %s %s %s\n',cluster.np*cluster.node,cluster.codepath,solution,[cluster.executionpath '/' dirname],modelname);
				 % NOTE: new version 
				 % not requires np processor and machine file. Execution of ISSM is operated under PBS. 
				 fprintf(fid,'mpiexec -genv MV2_ENABLE_AFFINITY 0 -np %d %s/issm.exe %s %s %s\n',...
					 cluster.np,cluster.codepath,solution,[cluster.executionpath '/' dirname],modelname);
			 else
				 % use machine file for obtaining cluster nodes
				 fprintf(fid,'mpiexec -machinefile %s -np %i %s/issm.exe %s %s %s\n',[cluster.executionpath '/' dirname '/simba.host'],cluster.np*cluster.node,cluster.codepath,solution,[cluster.executionpath '/' dirname],modelname);
			 end 
			 if ~io_gather, %concatenate the output files:
				 fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
			 end
			 fclose(fid);

			 %% generate machinefile "simba.host" file for SIMBA
			 %fid = fopen('simba.host','w');
			 %nodeorder = [1:cluster.numnodes, 0];

			 %% generate host of simba 
			 %for i = 1:cluster.numnodes-1,
			 %   fprintf(fid,'simba%02d:%2d\n',nodeorder(i),cluster.cpuspernode);
			 %end
			 %fclose(fid);

			 %in interactive mode, create a run file, and errlog and outlog file
			 if cluster.interactive,
				 fid=fopen([modelname '.run'],'w');
				 fprintf(fid,'mpiexec -np %i %s/issm.exe %s %s %s\n',cluster.np,cluster.codepath,solution,[cluster.executionpath '/' dirname],modelname);
				 if ~io_gather, %concatenate the output files:
					 fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
				 end
				 fclose(fid);
				 fid=fopen([modelname '.errlog'],'w');
				 fclose(fid);
				 fid=fopen([modelname '.outlog'],'w');
				 fclose(fid);
			 end
		 end %}}}
		 function UploadQueueJob(cluster,modelname,dirname,filelist)% {{{
			 %compress the files into one zip.
			 compressstring=['tar -zcf ' dirname '.tar.gz '];
			 for i=1:numel(filelist),
				 compressstring = [compressstring ' ' filelist{i}];
			 end
			 if cluster.interactive,
				 compressstring = [compressstring ' ' modelname '.errlog ' modelname '.outlog '];
			 end

			 % add machine file for SIMBA operation
			 %compressstring = [compressstring ' ' 'simba.host'];
			 system(compressstring);

			 if cluster.verbose
				 disp('hpc_simba2: uploading input file and queueing script');
				 if exist(sprintf('%s.tar.gz',dirname))
					fprintf('hpc_simba2: -- we find %s.tar.gz\n',dirname);
				 end
				 fprintf('hpc_simba2: -- remote hostname: %s\n',cluster.name);
			 end
			 
			 % use rclone for upload model.
			 isrclone = 0;
			 if isrclone,
				 command = ['ssh simba20 "rclone copy simba00:' pwd '/' dirname '.tar.gz ' cluster.executionpath '"'];
				 system(command);
			 else
				 if strfind(cluster.executionpath,'/data2/')
					 launchcommand = ['cp ' pwd '/' dirname '.tar.gz ' cluster.executionpath];
					 issmssh(oshostname,cluster.login,cluster.port,launchcommand);
				 else % using scpout for exporting data...
					 issmscpout(cluster.name,cluster.executionpath,cluster.login,cluster.port,{[dirname '.tar.gz']});
				 end
			 end
		 end %}}}
		 function LaunchQueueJob(cluster,modelname,dirname,filelist,restart,batch)% {{{
			 if ~isempty(restart)
				 if cluster.isqsub 
					 launchcommand=['cd ' cluster.executionpath ' && cd ' dirname ' && hostname && qsub ' modelname '.queue '];
				 else
					 launchcommand=['cd ' cluster.executionpath ' && cd ' dirname ' && hostname && source' modelname '.queue '];
				 end
			 else
				 if cluster.isqsub 
					 launchcommand=['cd ' cluster.executionpath ' && rm -rf ./' dirname ' && mkdir ' dirname ...
						 ' && cd ' dirname ' && mv ../' dirname '.tar.gz ./ && tar -zxf ' dirname '.tar.gz  && hostname && qsub ' modelname '.queue '];
				 else
					 launchcommand=['cd ' cluster.executionpath ' && rm -rf ./' dirname ' && mkdir ' dirname ...
						 ' && cd ' dirname ' && mv ../' dirname '.tar.gz ./ && tar -zxf ' dirname '.tar.gz  && hostname && source ' modelname '.queue '];
				 end
			 end
			 
			 if cluster.verbose
				 fprintf('check simulation at\n');
				 fprintf('%s/%s\n',cluster.executionpath,dirname);
			 end

			 % launch simulation
			 issmssh(cluster.name,cluster.login,cluster.port,launchcommand);
		 end %}}}
		 function Download(cluster,dirname,filelist)% {{{
			 %copy files from cluster to current directory
			 directory=[cluster.executionpath '/' dirname '/'];
			
			 isrclone = 0;
			 if isrclone,
				 command = ['ssh simba20 "rclone copy ' directory ' simba00:' pwd '/'];
				 for i = 1:length(filelist)
					 command = [command ' --include ' filelist{i}];
				 end
				 if cluster.verbose,
					 command = [command ' --progress '];
				 end
				 command = [command '"'];
				 %assignin('base','command',command); 

				 disp('download outputs....');
				 system(command);
				 disp('download outputs.... done...');
			 else
				 % In case of "SIMBA" machine, it shares the data storage. If following directory is set we do not need to download load file with using "issmscpin". Directly read file...
				 if 1==strfind(cluster.executionpath,'/data2/')
					 % tricky part for assign "cluster.name" as "simba00".
					 disp('hpc_simba2: download file from /data2/ directory!!! Not use scpin.');
					 issmscpin('simba00',cluster.login,cluster.port,directory,filelist);
				 else
					 issmscpin(cluster.name,cluster.login,cluster.port,directory,filelist);
				 end
			 end
		 end %}}}
	end
end
