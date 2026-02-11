%FRONTERA cluster class definition
%
%   Usage:
%      cluster=frontera();
%      cluster=frontera('np',3);
%      cluster=frontera('np',3,'login','username');

classdef frontera
	properties (SetAccess=public)
		% {{{
		name          = 'frontera'
		login         = '';
		modules        = {'intel/23.1.0', 'impi/21.9.0', 'petsc/3.21'};
		numnodes      = 4;
		cpuspernode   = 56;
		port          = 0;
		queue         = 'normal';
		codepath      = '';
		executionpath = '';
		interactive   = 0;
		time          = 5;
		email         = '';
	end
	%}}}
	methods
		function cluster=frontera(varargin) % {{{

			%initialize cluster using default settings if provided
			if (exist('frontera_settings')==2), frontera_settings; end

			%use provided options to change fields
			cluster=AssignObjectFields(pairoptions(varargin{:}),cluster);

		end
		%}}}
		function disp(cluster) % {{{
			%  display the object
			disp(sprintf('class ''%s'' object ''%s'' = ',class(cluster),inputname(1)));
			disp(sprintf('    name: %s',cluster.name));
			disp(sprintf('    login: %s',cluster.login));
			disp(sprintf('    modules: %s',strjoin(cluster.modules,', ')));
			disp(sprintf('    port: %i',cluster.port));
			disp(sprintf('    numnodes: %i',cluster.numnodes));
			disp(sprintf('    cpuspernode: %i',cluster.cpuspernode));
			disp(sprintf('    np: %i',cluster.nprocs()));
			disp(sprintf('    port: %i',cluster.port));
			disp(sprintf('    queue: %s',cluster.queue));
			disp(sprintf('    codepath: %s',cluster.codepath));
			disp(sprintf('    executionpath: %s',cluster.executionpath));
			disp(sprintf('    interactive: %i',cluster.interactive));
			disp(sprintf('    time: %i',cluster.time));
			disp(sprintf('    email: %s',cluster.email));
		end
		%}}}
		function numprocs=nprocs(cluster) % {{{
			%compute number of processors
			numprocs=cluster.numnodes*cluster.cpuspernode;
		end
		%}}}
		function md = checkconsistency(cluster,md,solution,analyses) % {{{

			%https://frontera-portal.tacc.utexas.edu/user-guide/running/
			available_queues={'flex','normal','development', 'small'};
			queue_requirements_time=[48 48 2 48];
			queue_requirements_np=[128*56 512*56 40*56 2*56];

			QueueRequirements(available_queues,queue_requirements_time,queue_requirements_np,cluster.queue,cluster.nprocs(),cluster.time)

			%Miscellaneous
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
			fprintf(fid,'mpiexec -np %i %s/kriging.exe %s %s\n',cluster.nprocs(),cluster.codepath,[cluster.executionpath '/' modelname],modelname);
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

			%write queuing script 
			fid=fopen([modelname '.queue'],'w');

			fprintf(fid,'#!/bin/bash\n');
			fprintf(fid,'#SBATCH -J %s \n',modelname);
			fprintf(fid,'#SBATCH -p %s \n',cluster.queue);
			fprintf(fid,'#SBATCH -o %s.outlog \n',modelname);
			fprintf(fid,'#SBATCH -e %s.errlog \n',modelname);
			fprintf(fid,'#SBATCH -n %i \n',cluster.numnodes*max(cluster.nprocs()/cluster.numnodes,56));
			fprintf(fid,'#SBATCH -N %i \n',cluster.numnodes);
			fprintf(fid,'#SBATCH -t %02i:00:00 \n\n',ceil(cluster.time));
			if length(find(cluster.email=='@'))>0
				fprintf(fid,'#SBATCH --mail-user=%s \n',cluster.email);
				fprintf(fid,'#SBATCH --mail-type=all \n\n');
			end
			for i=1:numel(cluster.modules),
				fprintf(fid,['module load ' cluster.modules{i} '\n']);
			end

			if isdakota,
				fprintf(fid,'export KMP_AFFINITY="granularity=fine,compact,verbose" \n\n');
			end


			fprintf(fid,'export PATH="$PATH:."\n\n');
			fprintf(fid,'export ISSM_DIR="%s/../"\n',cluster.codepath); %FIXME
			fprintf(fid,'source $ISSM_DIR/etc/environment.sh\n');       %FIXME
			fprintf(fid,'cd %s/%s\n\n',cluster.executionpath,dirname);
			fprintf(fid,'ibrun -np %i %s/%s %s %s %s\n',cluster.nprocs(),cluster.codepath,executable,solution,[cluster.executionpath '/' dirname],modelname);
			if ~io_gather, %concatenate the output files:
				fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
			end

			fclose(fid);

			%in interactive mode, create a run file, and errlog and outlog file
			if cluster.interactive,
				fid=fopen([modelname '.run'],'w');
				fprintf(fid,'ibrun -np %i %s/%s %s %s %s\n',cluster.nprocs(),executable,cluster.codepath,solution,[cluster.executionpath '/' dirname],modelname);
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
			if cluster.interactive,
				compressstring = [compressstring ' ' modelname '.errlog ' modelname '.outlog '];
			end
			system(compressstring);

			%upload input files
			issmscpout(cluster.name,cluster.executionpath,cluster.login,cluster.port,{[dirname '.tar.gz']});

		end %}}}
		function LaunchQueueJob(cluster,modelname,dirname,filelist,restart,batch) % {{{

			%Execute Queue job
			if ~isempty(restart)
				launchcommand=['cd ' cluster.executionpath ' && cd ' dirname ' && hostname && sbatch ' modelname '.queue '];
			else
				launchcommand=['cd ' cluster.executionpath ' && rm -rf ./' dirname ' && mkdir ' dirname ...
					' && cd ' dirname ' && mv ../' dirname '.tar.gz ./ && tar -zxf ' dirname '.tar.gz  && hostname && sbatch ' modelname '.queue '];
			end
			issmssh(cluster.name,cluster.login,cluster.port,launchcommand);
		end %}}}
		function Download(cluster,dirname,filelist) % {{{

			%copy files from cluster to current directory
			directory=[cluster.executionpath '/' dirname '/'];
			issmscpin(cluster.name,cluster.login,cluster.port,directory,filelist);

		end %}}}
	end
end
