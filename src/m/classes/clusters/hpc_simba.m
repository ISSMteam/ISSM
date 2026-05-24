%HPC_SIMBA cluster class definition
%
%   Usage:
%      cluster=hpc();
%      cluster=hpc('np',3);
%      cluster=hpc('np',3,'login','username');

classdef hpc_simba
	properties (SetAccess=public)
		% {{{
		name='simba00'
		login='inwoo';
		numnodes=1;
		cpuspernode=36; 
		np=4;
		port=7777;
		queue='pub64';
		codepath='/home/DATA/externalpackagesissm/issm/trunk-jpl/bin';
		executionpath='/home/DATA/externalpackages/issm/trunk-jpl/execution/';
		interactive=0;
	end
	%}}}
	methods
		function cluster=hpc_simba(varargin) % {{{

			%initialize cluster using default settings if provided
			if (exist('hpc_settings')==2), hpc_settings; end

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
			disp(sprintf('    np: %i',cluster.np));
			disp(sprintf('    port: %i',cluster.port));
			disp(sprintf('    queue: %s',cluster.queue));
			disp(sprintf('    codepath: %s',cluster.codepath));
			disp(sprintf('    executionpath: %s',cluster.executionpath));
			disp(sprintf('    interactive: %i',cluster.interactive));
		end
		%}}}
		function md = checkconsistency(cluster,md,solution,analyses) % {{{

			available_queues={'pub64','free64','free48','free*,pub64','free*'};
			queue_requirements_time=[Inf Inf Inf Inf Inf];
			queue_requirements_np=[64 64 48 48 48];

			QueueRequirements(available_queues,queue_requirements_time,queue_requirements_np,cluster.queue,cluster.np,1)

			%Miscellaneous
			if isempty(cluster.login), md = checkmessage(md,'login empty'); end
			if isempty(cluster.codepath), md = checkmessage(md,'codepath empty'); end
			if isempty(cluster.executionpath), md = checkmessage(md,'executionpath empty'); end

		end
		%}}}
		function BuildKrigingQueueScript(cluster, md, filename) % {{{

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
			if ~io_gather %concatenate the output files:
				fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
			end
			fclose(fid);
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
			if 0
				fprintf(fid,'#!/bin/bash\n');
				fprintf(fid,'#$ -N %s\n',modelname);
				fprintf(fid,'#$ -q %s \n',cluster.queue);
				fprintf(fid,'#$ -pe one-node-mpi 2-64\n');
				fprintf(fid,'#$ -R y\n');
				%fprintf(fid,'#$ -m beas\n');
				fprintf(fid,'#$ -o %s.outlog \n',modelname);
				fprintf(fid,'#$ -e %s.errlog \n\n',modelname);
			else % FIXME : This script for SIMBA-Kopri 
				fprintf(fid,'#!/bin/bash\n');
				%fprintf(fid,'#PBS -q workq\n');
				fprintf(fid,'#PBS -N %s\n',modelname);
				%fprintf(fid,'#PBS -l nodes=simba01:ppn=%d\n',...
				%	ceil(cluster.np/cluster.cpuspernode), cluster.np);
				fprintf(fid,'#PBS -l nodes=simba01:ppn=%d\n', cluster.np);
				fprintf(fid,'#PBS -o %s.outlog \n',modelname);
				fprintf(fid,'#PBS -e %s.errlog \n\n',modelname);
				fprintf(fid,'\n');
			end
			fprintf(fid,'export ISSM_DIR="%s/../"\n',cluster.codepath); %FIXME
			fprintf(fid,'source $ISSM_DIR/etc/environment.sh\n');       %FIXME
			fprintf(fid,'\n');
			fprintf(fid,'module load intel18/impi-18\n');
			fprintf(fid,'\n');
			fprintf(fid,'cd %s/%s\n\n',cluster.executionpath,dirname);
			fprintf(fid,'mpirun -np %i %s/%s %s %s %s\n',cluster.np,cluster.codepath,executable,solution,[cluster.executionpath '/' dirname],modelname);
			if ~io_gather, %concatenate the output files:
				fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
			end
			fclose(fid);

			%in interactive mode, create a run file, and errlog and outlog file
			if cluster.interactive
				fid=fopen([modelname '.run'],'w');
				fprintf(fid,'mpiexec -np %i %s/%s %s %s %s\n',cluster.np,cluster.codepath,executable,solution,[cluster.executionpath '/' dirname],modelname);
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
			if ~isempty(restart)
				launchcommand=['cd ' cluster.executionpath ' && cd ' dirname ' && hostname && qsub ' modelname '.queue '];
			else
				launchcommand=['cd ' cluster.executionpath ' && rm -rf ./' dirname ' && mkdir ' dirname ...
					' && cd ' dirname ' && mv ../' dirname '.tar.gz ./ && tar -zxf ' dirname '.tar.gz  && hostname && qsub ' modelname '.queue '];
			end
			issmssh(cluster.name,cluster.login,cluster.port,launchcommand);
		end %}}}
		function Download(cluster,dirname,filelist) % {{{
			cluster_defaults.Download(cluster,dirname,filelist);
		end %}}}
	end
end
