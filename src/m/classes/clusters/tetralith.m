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
		function BuildKrigingQueueScript(cluster,modelname,solution,io_gather,isvalgrind,isgprof) % {{{

			if(isvalgrind), disp('valgrind not supported by cluster, ignoring...'); end
			if(isgprof),    disp('gprof not supported by cluster, ignoring...'); end

			%compute number of processors
% 			cluster.np=cluster.numnodes*cluster.cpuspernode;
%			nprocs(cluster);%=cluster.numnodes*cluster.cpuspernode;

			%write queuing script 
			fid=fopen([modelname '.queue'],'w');
			fprintf(fid,'#!/bin/bash\n');
			fprintf(fid,'#\n');
			fprintf(fid,'#SBATCH --job-name=%s\n',modelname);
% 			fprintf(fid,'#SBATCH -p %s \n',cluster.partition);
			fprintf(fid,'#SBATCH -A %s \n',cluster.accountname);
% 			fprintf(fid,'#SBATCH --mail-type=ALL\n');
			fprintf(fid,'#SBATCH -N %i -n %i\n',cluster.numnodes,cluster.cpuspernode);
			%calculate walltime in hh:mm:ss format
			walltime=datestr(cluster.time/(60*24),'HH:MM:SS')
			fprintf(fid,'#SBATCH -t %s\n',walltime); %walltime should be in hh:mm:ss
			fprintf(fid,'#SBATCH --mem=%i\n',cluster.mem);
			fprintf(fid,'#SBATCH -o %s.outlog \n',modelname);
			fprintf(fid,'#SBATCH -e %s.errlog \n\n',modelname);
% 			fprintf(fid,'module load intelcomp/17.0.0\n') %module load not recommended within job script at Tetralith
% 			fprintf(fid,'module load mpt/2.14\n')
% 			fprintf(fid,'module load petsc/3.7.4d\n')
% 			fprintf(fid,'module load parmetis/4.0.3\n') 
% 			fprintf(fid,'module load mumps/5.0.2\n')
% 			fprintf(fid,'export ISSM_DIR="%s"\n',cluster.codepath); %FIXME
			fprintf(fid,'export ISSM_DIR="%s/../"\n',cluster.codepath); %FIXME
			fprintf(fid,'source $ISSM_DIR/etc/environment.sh\n');       %FIXME
			fprintf(fid,'cd %s/%s\n\n',cluster.executionpath,dirname);
% 			fprintf(fid,'mpiexec -np %i %s/%s %s %s %s\n',cluster.nprocs(),cluster.codepath,executable,solution,[cluster.executionpath '/' dirname],modelname);
% 			fprintf(fid,'mpiexec_mpt -np %i %s/%s %s %s %s\n',cluster.np,cluster.codepath,executable,solution,[cluster.executionpath '/' dirname],modelname);
			fprintf(fid,'mpiexec -np %i %s/issm.exe %s %s %s\n',cluster.nprocs(),cluster.codepath,solution,[cluster.executionpath '/' dirname],modelname);
% 			fprintf(fid,'mpirun -np %i %s/issm.exe %s %s %s\n',cluster.nprocs(),cluster.codepath,solution,[cluster.executionpath '/' dirname],modelname);
			if ~io_gather, %concatenate the output files:
				fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
			end
			fclose(fid);
		end
		%}}}
		function BuildQueueScript(cluster,dirname,modelname,solution,io_gather,isvalgrind,isgprof,isdakota,isoceancoupling) % {{{

			if(isvalgrind), disp('valgrind not supported by cluster, ignoring...'); end
			if(isgprof),    disp('gprof not supported by cluster, ignoring...'); end

			executable='issm.exe';
			if isdakota,
				version=IssmConfig('_DAKOTA_VERSION_'); version=str2num(version(1:3));
				if (version>=6),
					executable='issm_dakota.exe';
				end
			end
			if isoceancoupling,
				executable='issm_ocean.exe';
			end

			%compute number of processors
% 			cluster.np=cluster.numnodes*cluster.cpuspernode;
			nprocs(cluster);%=cluster.numnodes*cluster.cpuspernode;
% 			shortname = substring(modelname,1,min(12,length(modelname)));

			%write queuing script 
			fid=fopen([modelname '.queue'],'w');
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
			fprintf(fid,'mpiexec -np %i %s/issm.exe %s %s %s\n',cluster.nprocs(),cluster.codepath,solution,[cluster.executionpath '/' dirname],modelname);
% 			fprintf(fid,'mpirun -np %i %s/issm.exe %s %s %s\n',cluster.nprocs(),cluster.codepath,solution,[cluster.executionpath '/' dirname],modelname);
% 			fprintf(fid,'mpiexec_mpt -np %i %s/%s %s %s %s\n',cluster.nprocs(),cluster.codepath,executable,solution,[cluster.executionpath '/' dirname],modelname);

			if ~io_gather, %concatenate the output files:
				fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
			end
			fclose(fid);

			%in interactive mode, create a run file, and errlog and outlog file
			if cluster.interactive,
				fid=fopen([modelname '.run'],'w');
				fprintf(fid,'mpiexec_mpt -np %i %s/issm.exe %s %s %s\n',cluster.nprocs(),cluster.codepath,solution,[cluster.executionpath '/' dirname],modelname);
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
		function UploadQueueJob(cluster,modelname,dirname,filelist) % {{{

			%compress the files into one zip.
			compressstring=['tar -zcf ' dirname '.tar.gz '];
			for i=1:numel(filelist),
				compressstring = [compressstring ' ' filelist{i}];
			end
			system(compressstring);
			%upload input files
			directory=cluster.executionpath;
% 			issmbbftpout(cluster.name,directory,cluster.login,cluster.port,cluster.numstreams,{[dirname '.tar.gz']});
			issmscpout(cluster.name,directory,cluster.login,cluster.port,{[dirname '.tar.gz']});

		end
		%}}}
		function LaunchQueueJob(cluster,modelname,dirname,filelist,restart,batch) % {{{

			%Execute Queue job
			launchcommand=['cd ' cluster.executionpath ' && rm -rf ./' dirname ' && mkdir ' dirname ...
				' && cd ' dirname ' && mv ../' dirname '.tar.gz ./ && tar -zxf ' dirname '.tar.gz  && hostname && sbatch ' modelname '.queue '];

			issmssh(cluster.name,cluster.login,cluster.port,launchcommand);
		end %}}}
		function Download(cluster,dirname,filelist) % {{{

			%copy files from cluster to current directory
			directory=[cluster.executionpath '/' dirname '/'];
			issmscpin(cluster.name,cluster.login,cluster.port,directory,filelist);

		end %}}}
	end
end
