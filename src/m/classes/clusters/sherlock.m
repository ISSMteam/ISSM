%SHERLOCK cluster class definition
%
%   Usage:
%      cluster=sherlock();
%      cluster=sherlock('np',3);
%      cluster=sherlock('np',3,'login','username');

classdef sherlock
	properties (SetAccess=public)
		% {{{
		name          = 'sherlock'
		login         = '';
		numnodes      = 1;
		cpuspernode   = 24;
		port          = 0;
		codepath      = '';
		executionpath = '';
		interactive   = 0;
		time          = 30;
		memory        = 2;
	end
	%}}}
	methods
		function cluster=sherlock(varargin) % {{{

			%initialize cluster using default settings if provided
			if (exist('sherlock_settings')==2), sherlock_settings; end

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

			%Miscellaneous
			if isempty(cluster.login), md = checkmessage(md,'login empty'); end
			if isempty(cluster.codepath), md = checkmessage(md,'codepath empty'); end
			if isempty(cluster.executionpath), md = checkmessage(md,'executionpath empty'); end

		end
		%}}}
		function BuildQueueScript(cluster, md, filename, executable) % {{{

			%Get variables from md
			dirname   = md.private.runtimename;
			modelname = md.miscellaneous.name;
			solution  = md.private.solution;
			io_gather = md.settings.io_gather;

			%checks
			if(md.debug.valgrind) disp('valgrind not supported by cluster, ignoring...'); end
			if(md.debug.gprof)    disp('gprof not supported by cluster, ignoring...'); end

			%write queuing script
			fid=fopen(filename, 'w');
			fprintf(fid,'#!/bin/bash\n');
			fprintf(fid,'#SBATCH --job-name=%s\n',modelname);
			fprintf(fid,'#SBATCH -N %i -n %i\n',cluster.numnodes,cluster.cpuspernode);
			fprintf(fid,'#SBATCH --time=%i\n',cluster.time*60); %walltime is in seconds.
			fprintf(fid,'#SBATCH --mem-per-cpu=%igb\n',cluster.memory);
			fprintf(fid,'#SBATCH -o %s.outlog \n',modelname);
			fprintf(fid,'#SBATCH -e %s.errlog \n\n',modelname);
			fprintf(fid,'export ISSM_DIR="%s/../"\n',cluster.codepath); %FIXME
			fprintf(fid,'source $ISSM_DIR/etc/environment.sh\n');       %FIXME
			fprintf(fid,'cd %s/%s\n\n',cluster.executionpath,dirname);
			fprintf(fid,'mpiexec -n %i %s/%s %s %s %s\n',cluster.nprocs(),cluster.codepath,executable,solution,[cluster.executionpath '/' dirname],modelname);
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
