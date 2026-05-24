%GREENPLANET cluster class definition
%
%   Usage:
%      cluster=greenplanet();
%      cluster=greenplanet('np',3);
%      cluster=greenplanet('np',3,'login','username');

classdef greenplanet
	properties (SetAccess=public)
		% {{{
		name          = 'greenplanet'
		login         = '';
		numnodes      = 20;
		cpuspernode   = 8;
		port          = 8000;
		queue         = 'c6145';
		codepath      = '';
		executionpath = '';
		interactive   = 0;
		time          = 24*60;
		memory        = 2;
	end
	%}}}
	methods
		function cluster=greenplanet(varargin) % {{{

			%initialize cluster using default settings if provided
			if (exist('greenplanet_settings')==2), greenplanet_settings; end

			%use provided options to change fields
			cluster=AssignObjectFields(pairoptions(varargin{:}),cluster);
		end
		%}}}
		function disp(cluster) % {{{
			%  display the object
			disp(sprintf('class ''%s'' object ''%s'' = ',class(cluster),inputname(1)));
			disp(sprintf('    name: %s',cluster.name));
			disp(sprintf('    login: %s',cluster.login));
			disp(sprintf('    numnodes: %i',cluster.numnodes));
			disp(sprintf('    cpuspernode: %i',cluster.cpuspernode));
			disp(sprintf('    np: %i',cluster.nprocs()));
			disp(sprintf('    port: %i',cluster.port));
			disp(sprintf('    queue: %s',cluster.queue));
			disp(sprintf('    codepath: %s',cluster.codepath));
			disp(sprintf('    executionpath: %s',cluster.executionpath));
			disp(sprintf('    interactive: %i',cluster.interactive));
			disp(sprintf('    time: %i',cluster.time));
			disp(sprintf('    memory: %i',cluster.memory));
		end
		%}}}
		function numprocs=nprocs(cluster) % {{{
			%compute number of processors
			numprocs=cluster.numnodes*cluster.cpuspernode;
		end
		%}}}
		function md = checkconsistency(cluster,md,solution,analyses) % {{{

			available_queues={'brd2.4','has2.5','ilg2.3','m-c1.9','m-c2.2','nes2.8','sib2.9','m2090','default'};
			queue_requirements_time=[Inf Inf Inf Inf Inf Inf Inf Inf Inf];
			queue_requirements_np=[80 80 80 80 80 80 80 80 80];

			QueueRequirements(available_queues,queue_requirements_time,queue_requirements_np,cluster.queue,cluster.nprocs(),1)

			%Miscellaneous
			if isempty(cluster.login), md = checkmessage(md,'login empty'); end
			if isempty(cluster.codepath), md = checkmessage(md,'codepath empty'); end
			if isempty(cluster.executionpath), md = checkmessage(md,'executionpath empty'); end

		end
		%}}}
		function BuildQueueScript(cluster, md, filename, executable) % {{{

         %Get variables from md
         dirname         = md.private.runtimename;
         modelname       = md.miscellaneous.name;
         solution        = md.private.solution;
         io_gather       = md.settings.io_gather;

         %checks
			if(md.debug.valgrind) disp('valgrind not supported by cluster, ignoring...'); end
			if(md.debug.gprof)    disp('gprof not supported by cluster, ignoring...'); end

			%write queuing script 
			fid=fopen(filename, 'w');
			fprintf(fid,'#!/bin/bash\n');
			fprintf(fid,'#SBATCH --job-name=%s\n',modelname);
			fprintf(fid,'#SBATCH --partition=%s',cluster.queue{1});
			for i=2:length(cluster.queue) 
				fprintf(fid,',%s',cluster.queue{i});
			end
			fprintf(fid,'\n');
			%fprintf(fid,'#SBATCH -N %i -n %i\n',cluster.numnodes,cluster.cpuspernode);
			%fprintf(fid,'#SBATCH --mem-per-cpu=%igb\n',cluster.memory);
			fprintf(fid,'#SBATCH --nodes=%i\n',cluster.numnodes); % in general, just 1
			fprintf(fid,'#SBATCH --ntasks=%i\n',cluster.cpuspernode); % in general, just 1
			fprintf(fid,'#SBATCH --cpus-per-task=%i\n',1); 
			fprintf(fid,'#SBATCH --mem=%igb\n',cluster.memory); % minimum total node memory required
			fprintf(fid,'#SBATCH --time=%s\n',datestr(cluster.time/24,'HH:MM:SS')); %walltime is in HH:MM:SS format. cluster.time is in hour
			fprintf(fid,'#SBATCH -o %s.outlog \n',modelname);
			fprintf(fid,'#SBATCH -e %s.errlog \n\n',modelname);
			fprintf(fid,'export ISSM_DIR="%s/../"\n',cluster.codepath); %FIXME
			fprintf(fid,'source $ISSM_DIR/etc/environment.sh\n');       %FIXME
			fprintf(fid,'cd %s/%s\n\n',cluster.executionpath,dirname);
			fprintf(fid,'mpiexec -n %i %s/%s %s %s %s\n',cluster.nprocs(),cluster.codepath,executable,solution,[cluster.executionpath '/' dirname],modelname);
			if ~io_gather, %concatenate the output files:
				fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
			end
			fclose(fid);

			%in interactive mode, create a run file, and errlog and outlog file
			if cluster.interactive
				fid=fopen([modelname '.run'],'w');
				fprintf(fid,'mpiexec -n %i %s/%s %s %s %s\n',cluster.nprocs(),cluster.codepath,executable,solution,[cluster.executionpath '/' dirname],modelname);
				if ~io_gather, %concatenate the output files:
					fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
				end
				fclose(fid);
				fid=fopen([modelname '.errlog'],'w'); fclose(fid);
				fid=fopen([modelname '.outlog'],'w'); fclose(fid);
			end
		end %}}}
		function UploadQueueJob(cluster,modelname,dirname,filelist) % {{{
			cluster_defaults.UploadQueueJob(cluster,modelname,dirname,filelist);
		end %}}}
		function LaunchQueueJob(cluster,modelname,dirname,filelist,restart,batch) % {{{
			cluster_defaults.LaunchQueueJobSbatch(cluster,modelname,dirname,filelist,restart,batch, 2);
		end %}}}
		function Download(cluster,dirname,filelist) % {{{
			cluster_defaults.Download(cluster,dirname,filelist);
		end %}}}
	end
end
