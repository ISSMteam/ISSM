%AURORA cluster class definition
%
%   Usage:
%      cluster=aurora();
%      cluster=aurora('np',3);
%      cluster=aurora('np',3,'login','username');

classdef aurora
	properties (SetAccess=public)
		% {{{
		name='aurora.jpl.nasa.gov'
		login='username';
		numnodes = 1;
		cpuspernode = 24;
		modules = {'intel/cluster-toolkit-2013.5.192'};
		port=1070;
		queue='shortq';
		time=3*60;
		codepath='~/issm/trunk/';
		executionpath='~/issm/trunk/execution/';
		mpipath='/opt/intel/impi/4.1.3/intel64/bin/';
		%}}}
	end
	methods
		function cluster=aurora(varargin) % {{{

			%initialize cluster using default settings if provided
			if (exist('aurora_settings')==2), aurora_settings; end

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
			disp(sprintf('    modules: %s',strjoin(cluster.modules,', ')));
			disp(sprintf('    port: %i',cluster.port));
			disp(sprintf('    queue: %s',cluster.queue));
			disp(sprintf('    time: %i',cluster.time));
			disp(sprintf('    codepath: %s',cluster.codepath));
			disp(sprintf('    executionpath: %s',cluster.executionpath));
			disp(sprintf('    mpipath: %s',cluster.mpipath));
		end
		%}}}
		function numprocs=nprocs(cluster) % {{{
			%compute number of processors
			numprocs=cluster.numnodes*cluster.cpuspernode;
		end
		%}}}
		function md = checkconsistency(cluster,md,solution,analyses) % {{{

			available_queues={'debugq','shortq','mediumq','longq','verylongq'};
			queue_requirements_time=[60*1 60*3 60*12 60*48 60*192];
			queue_requirements_np=[16 256 256 128 128];

			QueueRequirements(available_queues,queue_requirements_time,queue_requirements_np,cluster.queue,cluster.numnodes.*cluster.cpuspernode,cluster.time)
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
			fprintf(fid,'#PBS -l select=%i:ncpus=%i\n',cluster.numnodes,cluster.cpuspernode);
			fprintf(fid,'#PBS -N %s\n',modelname);
			fprintf(fid,'#PBS -l walltime=%i\n',cluster.time*60); %walltime is in seconds.
			fprintf(fid,'#PBS -q %s\n',cluster.queue);
			fprintf(fid,'#PBS -o %s.outlog \n',modelname);
			fprintf(fid,'#PBS -e %s.errlog \n',modelname);
			fprintf(fid,'source /usr/share/Modules/init/bash\n');
			for i=1:numel(cluster.modules), fprintf(fid,['module load ' cluster.modules{i} '\n\n']); end
			fprintf(fid,'export PATH="$PATH:."\n\n');
			fprintf(fid,'export MPI_GROUP_MAX=64\n\n');
			fprintf(fid,'export MPI_UNBUFFERED_STDIO=true\n\n');
			fprintf(fid,'export ISSM_DIR="%s/../"\n',cluster.codepath);
			fprintf(fid,'source $ISSM_DIR/etc/environment.sh\n');      
			fprintf(fid,'export PATH="$PATH:.:%s"\n',cluster.mpipath);
			fprintf(fid,'export PBS_O_WORKDIR=%s\n',[cluster.executionpath '/' dirname]);
			fprintf(fid,'cd $PBS_O_WORKDIR\n');
			fprintf(fid,'mpirun -n %i %s/%s %s %s %s',cluster.nprocs(),cluster.codepath,executable,solution,[cluster.executionpath '/' dirname],modelname);
			fclose(fid);

		end
		%}}}
		function UploadQueueJob(cluster,modelname,dirname,filelist) % {{{

			%compress the files into one zip.
			compressstring=['tar -zcf ' dirname '.tar.gz '];
			for i=1:numel(filelist),
				compressstring = [compressstring ' ' filelist{i}];
			end
			system(compressstring);

			%upload input files
			issmscpout(cluster.name,cluster.executionpath,cluster.login,cluster.port,{[dirname '.tar.gz']});

		end %}}}
		function LaunchQueueJob(cluster,modelname,dirname,filelist,restart,batch) % {{{

			if ~isempty(restart)
				launchcommand=['cd ' cluster.executionpath ' && cd ' dirname ' && qsub ' modelname '.queue '];
			else
				launchcommand=['cd ' cluster.executionpath ' && rm -rf ./' dirname ' && mkdir ' dirname ...
					' && cd ' dirname ' && mv ../' dirname '.tar.gz ./ && tar -zxf ' dirname '.tar.gz  && qsub ' modelname '.queue '];
			end

			%Execute Queue job
			issmssh(cluster.name,cluster.login,cluster.port,launchcommand);

		end %}}}
		function Download(cluster,dirname,filelist) % {{{

			%copy files from cluster to current directory
			directory=[cluster.executionpath '/' dirname '/'];
			issmscpin(cluster.name,cluster.login,cluster.port,directory,filelist);

		end %}}}
	end
end
