%MAUI cluster class definition
%
%   Usage:
%      cluster=maui();
%      cluster=maui('np',3);
%      cluster=maui('np',3,'login','username');

classdef maui
	properties (SetAccess=public)
		% {{{
		name           = 'maui'
		login          = '';
		numnodes       = 1;
		cpuspernode    = 8;
		port           = 0;
		projectaccount = '';
		partition      = 'nesi_research';
		codepath       = '';
		executionpath  = '';
		interactive    = 0;
		time           = 24*60;
		memory         = 2;
	end
	%}}}
	methods
		function cluster=maui(varargin) % {{{

			%initialize cluster using default settings if provided
			if (exist('maui_settings')==2), maui_settings; end

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
			disp(sprintf('    projectaccount: %s',cluster.projectaccount));
			disp(sprintf('    partition: %s',cluster.partition));
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

			available_partitions={'nesi_research'}
			partition_requirements_time=[24*60];
			partition_requirements_np=[80];

			QueueRequirements(available_partitions,partition_requirements_time,partition_requirements_np,cluster.partition,cluster.nprocs(),1)

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
			fprintf(fid,'#SBATCH --account=%s \n',cluster.projectaccount);
			fprintf(fid,'#SBATCH --partition=%s \n',cluster.partition);
			fprintf(fid,'#SBATCH --ntasks=%i  \n',cluster.numnodes);
			fprintf(fid,'#SBATCH --cpus-per-task=%i\n',cluster.cpuspernode);
			fprintf(fid,'#SBATCH --time=%i\n',cluster.time); %walltime is in minutes
			fprintf(fid,'#SBATCH --mem-per-cpu=%igb\n',cluster.memory);
			fprintf(fid,'#SBATCH -o %s.outlog \n',modelname);
			fprintf(fid,'#SBATCH -e %s.errlog \n\n',modelname);
			fprintf(fid,'export ISSM_DIR="%s/../"\n',cluster.codepath); %FIXME
			fprintf(fid,'source $ISSM_DIR/etc/environment.sh\n');       %FIXME
			fprintf(fid,'module swap PrgEnv-cray PrgEnv-intel\n');
			fprintf(fid,'cd %s/%s\n\n',cluster.executionpath,dirname);
			fprintf(fid,'srun -n %i %s/%s %s %s %s\n',cluster.nprocs(),cluster.codepath,executable,solution,[cluster.executionpath '/' dirname],modelname);
			if ~io_gather %concatenate the output files:
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
				fid=fopen([modelname '.errlog'],'w'); fclose(fid); fid=fopen([modelname '.outlog'],'w');
				fclose(fid);
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
