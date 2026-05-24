%TERALITH cluster class definition
%
%   Usage:
%      cluster=tetralith();
%      cluster=tetralith('np',3);
%      cluster=tetralith('np',3,'login','username');

classdef tetralith
	properties (SetAccess=public)
		% {{{
		name           = 'tetralith';
		login          = '';
		numnodes       = 2;
		cpuspernode    = 32;
		mem            = 2000; %MB
		queue          = 'normal';
		time           = 2*60; %[minutes]
		codepath       = '';
		executionpath  = '';
		interactive    = 0;
		port           = [];
		accountname    = '';
		% }}}
	end
	methods
		function cluster=tetralith(varargin) % {{{

			%initialize cluster using default settings if provided
			if (exist('tetralith_settings')==2), tetralith_settings; end

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
			disp(sprintf('    memory: %i',cluster.mem));
			disp(sprintf('    queue: %s',cluster.queue));
			disp(sprintf('    time: %i',cluster.time));
			disp(sprintf('    codepath: %s',cluster.codepath));
			disp(sprintf('    executionpath: %s',cluster.executionpath));
			disp(sprintf('    interactive: %i',cluster.interactive));
			disp(sprintf('    port: %s',cluster.port));
			disp(sprintf('    accountname: %s',cluster.accountname));
		end
		%}}}
		function md = checkconsistency(cluster,md,solution,analyses) % {{{

			available_queues={'normal','devel'};
			queue_requirements_time=[7*24*60 60]; 
			queue_requirements_np=[1024 2*32];

			QueueRequirements(available_queues,queue_requirements_time,queue_requirements_np,cluster.queue,cluster.nprocs(),1)

			%Miscellaneous
			if isempty(cluster.login), md = checkmessage(md,'login empty'); end
			if isempty(cluster.accountname), md = checkmessage(md,'accountname empty'); end
			if isempty(cluster.codepath), md = checkmessage(md,'codepath empty'); end
			if isempty(cluster.executionpath), md = checkmessage(md,'executionpath empty'); end

		end
		%}}}
		function numprocs=nprocs(self) % {{{
			%compute number of processors
			numprocs=self.numnodes*self.cpuspernode;
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
			fprintf(fid,'#\n');
			fprintf(fid,'#SBATCH --job-name=%s\n',modelname);
% 			fprintf(fid,'#SBATCH -p %s \n',cluster.partition);
			fprintf(fid,'#SBATCH -A %s \n',cluster.accountname);
% 			fprintf(fid,'#SBATCH --mail-type=ALL\n');
			fprintf(fid,'#SBATCH -N %i -n %i\n',cluster.numnodes,cluster.cpuspernode);
			%calculate walltime in hh:mm:ss format
			walltime=datestr(cluster.time/(60*24),'HH:MM:SS');
			fprintf(fid,'#SBATCH -t %s\n',walltime); %walltime should be in hh:mm:ss
			fprintf(fid,'#SBATCH --mem=%i\n',cluster.mem);
			fprintf(fid,'#SBATCH -o %s.outlog \n',modelname);
			fprintf(fid,'#SBATCH -e %s.errlog \n\n',modelname);
% 			fprintf(fid,'module load intelcomp/17.0.0\n') %module load not recommended within job script at Tetralith
% 			fprintf(fid,'module load mpt/2.14\n')
% 			fprintf(fid,'module load petsc/3.7.4d\n')
% 			fprintf(fid,'module load parmetis/4.0.3\n') 
% 			fprintf(fid,'module load mumps/5.0.2\n')
			fprintf(fid,'export ISSM_DIR="%s/../"\n',cluster.codepath); %FIXME
% 			fprintf(fid,'export ISSM_DIR="%s"\n',cluster.codepath); %FIXME
			fprintf(fid,'source $ISSM_DIR/etc/environment.sh\n');       %FIXME
			fprintf(fid,'cd %s/%s\n\n',cluster.executionpath,dirname);
% 			fprintf(fid,'mpiexec -np %i %s/%s %s %s %s\n',cluster.nprocs(),cluster.codepath,executable,solution,[cluster.executionpath '/' dirname],modelname);
			fprintf(fid,'mpiexec -np %i %s/%s %s %s %s\n',cluster.nprocs(),cluster.codepath,executable,solution,[cluster.executionpath '/' dirname],modelname);
% 			fprintf(fid,'mpirun -np %i %s/issm.exe %s %s %s\n',cluster.nprocs(),cluster.codepath,solution,[cluster.executionpath '/' dirname],modelname);
% 			fprintf(fid,'mpiexec_mpt -np %i %s/%s %s %s %s\n',cluster.nprocs(),cluster.codepath,executable,solution,[cluster.executionpath '/' dirname],modelname);

			if ~io_gather, %concatenate the output files:
				fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
			end
			fclose(fid);

			%in interactive mode, create a run file, and errlog and outlog file
			if cluster.interactive
				fid=fopen([filename '.run'],'w');
				fprintf(fid,'mpiexec_mpt -np %i %s/%s %s %s %s\n',cluster.nprocs(),cluster.codepath,executable,solution,[cluster.executionpath '/' dirname],modelname);
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
