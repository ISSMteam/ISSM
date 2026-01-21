%YELLOWSTONE cluster class definition
%
%   Usage:
%      cluster=yellowstone();
%      cluster=yellowstone('np',3);
%      cluster=yellowstone('np',3,'login','username');

classdef yellowstone
	properties (SetAccess=public)
		% {{{
		name           = 'yellowstone.ucar.edu'
		login          = '';
		modules        = {'ncarenv/1.0' 'ncarbinlibs/1.1' 'perlmods/5.0' 'gmake/4.1' 'python/2.7.7' 'all-python-libs' 'git/2.3.0' 'intel/15.0.3' 'mkl/11.1.2' 'esmf' 'esmf-6.3.0rp1-defio-mpi-O' 'netcdf-mpi/4.3.3.1' 'pnetcdf/1.6.1' 'ncarcompilers/1.0' 'cmake/3.0.2' 'matlab/R2015b' 'issm'};
		numnodes       = 1;
		cpuspernode    = 8;
		port           = 0;
		queue          = 'premium';
		time           = 12*60;
		processor      = 'sandy';
		codepath       = '';
		executionpath  = '';
		grouplist     = 'P93300301';
	end
	%}}}
	methods
		function cluster=yellowstone(varargin) % {{{

			%initialize cluster using default settings if provided
			if (exist('yellowstone_settings')==2), yellowstone_settings; end

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
			disp(sprintf('    numnodes: %i',cluster.numnodes));
			disp(sprintf('    cpuspernode: %i',cluster.cpuspernode));
			disp(sprintf('    np: %i',cluster.nprocs()));
			disp(sprintf('    port: %i',cluster.port));
			disp(sprintf('    queue: %s',cluster.queue));
			disp(sprintf('    time: %i',cluster.time));
			disp(sprintf('    processor: %s',cluster.processor));
			disp(sprintf('    codepath: %s',cluster.codepath));
			disp(sprintf('    executionpath: %s',cluster.executionpath));
			disp(sprintf('    grouplist: %s',cluster.grouplist));
		end
		%}}}
		function numprocs=nprocs(cluster) % {{{
			%compute number of processors
			numprocs=cluster.numnodes*cluster.cpuspernode;
		end
		%}}}
		function md = checkconsistency(cluster,md,solution,analyses) % {{{

			available_queues={'premium','regular'};
			queue_requirements_time=[12*60 12*650];
			queue_requirements_np=[16384 16384];

			QueueRequirements(available_queues,queue_requirements_time,queue_requirements_np,cluster.queue,cluster.nprocs(),cluster.time)

			%Miscellaneous
			if isempty(cluster.login), md = checkmessage(md,'login empty'); end
			if isempty(cluster.codepath), md = checkmessage(md,'codepath empty'); end
			if isempty(cluster.executionpath), md = checkmessage(md,'executionpath empty'); end
			if isempty(cluster.grouplist), md = checkmessage(md,'grouplist empty'); end

		end
		%}}}
		function BuildQueueScript(cluster,dirname,modelname,solution,io_gather,isvalgrind,isgprof,isdakota,isoceancoupling) % {{{

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
			fprintf(fid,'#!/bin/tcsh\n');
			fprintf(fid,'#BSUB -P %s\n',cluster.grouplist);
			fprintf(fid,'#BSUB -W %i:%i\n',floor(cluster.time/60),cluster.time-floor(cluster.time/60)*60);
			fprintf(fid,'#BSUB -n %i\n',cluster.nprocs());
			fprintf(fid,'#BSUB -J %s\n',modelname);
			fprintf(fid,'#BSUB -o %s.outlog \n',[cluster.executionpath '/' dirname '/' modelname]);
			fprintf(fid,'#BSUB -e %s.errlog \n',[cluster.executionpath '/' dirname '/' modelname]);
			fprintf(fid,'#BSUB -q %s\n',cluster.queue);

			fprintf(fid,'module purge\n');
			for i=1:length(cluster.modules),
				fprintf(fid,'module load %s\n',cluster.modules{i});
			end

			fprintf(fid,'setenv OMP_STACKSIZE 256M\n');
			fprintf(fid,'setenv MP_LABELIO yes\n');
			fprintf(fid,'setenv MP_INFOLEVEL 2\n');
			fprintf(fid,'setenv MP_SHARED_MEMORY yes\n');
			fprintf(fid,'setenv MP_EUILIB us\n');
			fprintf(fid,'setenv MP_STDOUTMODE unordered\n');
			fprintf(fid,'setenv MP_RC_USE_LMC yes\n');
			fprintf(fid,'setenv MP_MPILIB mpich2\n');
			fprintf(fid,'setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/opt/ibmhpc/pecurrent/mpich2/intel/lib64/\n');
			fprintf(fid,'setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/ncar/opt/intel/psxe-2015_update3/impi/5.0.3.048/lib64/\n');

			fprintf(fid,'cd %s/%s/\n\n',cluster.executionpath,dirname);

			fprintf(fid,'mpirun.lsf %s/%s %s %s %s\n',cluster.codepath,executable,solution,[cluster.executionpath '/' dirname],modelname);
			fclose(fid);

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

			issmscpout(cluster.name,directory,cluster.login,cluster.port,{[dirname '.tar.gz']});

		end
		%}}}
		function LaunchQueueJob(cluster,modelname,dirname,filelist,restart,batch) % {{{

			%launch command, to be executed via ssh
			launchcommand=['cd ' cluster.executionpath ' && rm -rf ./' dirname ' && mkdir ' dirname ...
						' && cd ' dirname ' && mv ../' dirname '.tar.gz ./ && tar -zxf ' dirname '.tar.gz  && bsub < ' modelname '.queue '];

			%Execute Queue job
			issmssh(cluster.name,cluster.login,cluster.port,launchcommand);
		end
		%}}}
		function Download(cluster,dirname,filelist) % {{{

			%copy files from cluster to current directory
			directory=[cluster.executionpath '/' dirname '/'];
			issmscpin(cluster.name,cluster.login,cluster.port,directory,filelist);

		end %}}}
	end
end
