%GADI cluster class definition
%
%   Usage:
%      cluster=gadi();
%      cluster=gadi('np',3);
%      cluster=gadi('np',3,'login','username');

classdef gadi 
	properties (SetAccess=public)
		% {{{
		name				= oshostname()
		login				= '';
		moduleload     = {};
		moduleuse      = {};
		numnodes       = 1;
		cpuspernode		= 48; 
		memory			= 40 % e.g., '40GB'
		port           = 0; % typical SSH port
		queue          = 'normal';
		time           = 60; % total minutes of walltime, e.g. 60 => 1hour
		processor      = ''; % not usually needed for Gadi
		srcpath        = '';
		extpkgpath     = '';
		codepath       = '';
		executionpath  = '';
		project        = '';
		storage        = '';
		interactive    = 0;
		bbftp          = 0;
		numstreams     = 8;
		hyperthreading = 0;
	end
	%}}}
	methods
		function cluster=gadi(varargin) % {{{

			%initialize cluster using default settings if provided
			if (exist('gdai_settings')==2), gadi_settings; end
			gadi_settings

			%use provided options to change fields
			cluster=AssignObjectFields(pairoptions(varargin{:}),cluster);
		end
		%}}}
		function disp(cluster) % {{{
			%  display the object
			disp(sprintf('class ''%s'' object ''%s'' = ',class(cluster),inputname(1)));
			disp(sprintf('    name: %s',cluster.name));
			disp(sprintf('    login: %s',cluster.login));
			disp(sprintf('    moduleuse: %s',strjoin(cluster.moduleuse,', ')));
			disp(sprintf('    moduleload: %s',strjoin(cluster.moduleload,', ')));
			disp(sprintf('    numnodes: %i',cluster.numnodes));
			disp(sprintf('    cpuspernode: %i',cluster.cpuspernode));
			disp(sprintf('    np: %i',cluster.nprocs()));
			disp(sprintf('    port: %i',cluster.port));
			disp(sprintf('    queue: %s',cluster.queue));
			disp(sprintf('    time (min): %.f',cluster.time));
			disp(sprintf('    processor: %s',cluster.processor));
			disp(sprintf('    srcpath: %s',cluster.srcpath));
			disp(sprintf('    extpkgpath: %s',cluster.extpkgpath));
			disp(sprintf('    codepath: %s',cluster.codepath));
			disp(sprintf('    executionpath: %s',cluster.executionpath));
			disp(sprintf('    project: %s',cluster.project));
			disp(sprintf('    storage: %s',cluster.storage));
			disp(sprintf('    interactive: %i',cluster.interactive));
			disp(sprintf('    bbftp: %i',cluster.bbftp));
			disp(sprintf('    numstreams: %i',cluster.numstreams));
			disp(sprintf('    hyperthreading: %i',cluster.hyperthreading));
		end
		%}}}
		function numprocs=nprocs(cluster) % {{{
			%compute number of processors
			numprocs=cluster.numnodes*cluster.cpuspernode;
		end
		%}}}
		function md = checkconsistency(cluster,md,solution,analyses) % {{{

			queuestruct = struct('normal',[48*60, 3072],... % e.g. up to 48h, 3072 cores
				'express',[2*60, 960],... % e.g. up to 2h, 960 cores
				'hugemem',[48*60, 3072]);
			available_queues = fieldnames(queuestruct);
			queue_requirements_time = zeros(length(available_queues),1);
			queue_requirements_np   = zeros(length(available_queues),1);
			for i = 1:length(available_queues)
				field = available_queues{i};
				queue_requirements_time(i) = queuestruct.(field)(1);
				queue_requirements_np(i) = queuestruct.(field)(2);
			end

			QueueRequirements(available_queues,queue_requirements_time,queue_requirements_np,cluster.queue,cluster.nprocs(),1)

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

			% Convert self.time (minutes) to hh:mm:ss
			hours  = floor(cluster.time/60);
			minutes= mod(cluster.time,60);
			walltime_str = sprintf('%02i:%02i:00',hours,minutes);

			%write queuing script 
			fid=fopen(filename, 'w');
			fprintf(fid,'#!/bin/bash\n');
			fprintf(fid,'#PBS -P %s\n',cluster.project);
			fprintf(fid,'#PBS -q %s\n',cluster.queue);
			fprintf(fid,'#PBS -l ncpus=%i\n',cluster.nprocs());
			fprintf(fid,'#PBS -l mem=%iGB\n',cluster.memory);
			fprintf(fid,'#PBS -l walltime=%s\n',walltime_str);
			fprintf(fid,'#PBS -l wd\n');
			fprintf(fid,'#PBS -j oe\n');
			fprintf(fid,'#PBS -l storage=%s\n',cluster.storage);
			fprintf(fid,'#PBS -m bea\n');
			fprintf(fid,'#PBS -o %s.outlog \n',modelname);
			fprintf(fid,'#PBS -e %s.errlog \n\n',modelname);

			fprintf(fid,'\n# Load modules as needed:\n');
			fprintf(fid,'module purge\n');
			for i=1:length(cluster.moduleuse);
				fprintf(fid,'module use %s\n',cluster.moduleuse{i});
			end
			for i=1:length(cluster.moduleload);
				fprintf(fid,'module load %s\n',cluster.moduleload{i});
			end

			fprintf(fid,'export ISSM_DIR="%s/../"\n',cluster.codepath); %FIXME
			fprintf(fid,'source $ISSM_DIR/etc/environment.sh\n');       %FIXME

			fprintf(fid,'\n# Switch to run directory (if not using -l wd):\n');
			fprintf(fid,'cd %s/%s\n\n',cluster.executionpath,dirname);

			fprintf(fid,'mpiexec -np %i %s/kriging.exe %s %s\n',cluster.nprocs(),cluster.codepath,[cluster.executionpath '/' modelname],modelname);
			if ~io_gather, %concatenate the output files:
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

			% Convert self.time (minutes) to hh:mm:ss
			hours  = floor(cluster.time/60);
			minutes= mod(cluster.time,60);
			walltime_str = sprintf('%02i:%02i:00',hours,minutes);

			%write queuing script 
			fid=fopen(filename, 'w');
			fprintf(fid,'#!/bin/bash\n');
			fprintf(fid,'#PBS -P %s\n',cluster.project);
			fprintf(fid,'#PBS -q %s\n',cluster.queue);
			fprintf(fid,'#PBS -l ncpus=%i\n',cluster.nprocs());
			fprintf(fid,'#PBS -l mem=%iGB\n',cluster.memory);
			fprintf(fid,'#PBS -l walltime=%s\n',walltime_str);
			fprintf(fid,'#PBS -l wd\n');
			fprintf(fid,'#PBS -j oe\n');
			fprintf(fid,'#PBS -l storage=%s\n',cluster.storage);
			fprintf(fid,'#PBS -m bea\n');
			fprintf(fid,'#PBS -o %s.outlog \n',modelname);
			fprintf(fid,'#PBS -e %s.errlog \n\n',modelname);

			fprintf(fid,'\n# Load modules as needed:\n');
			fprintf(fid,'module purge\n');
			for i=1:length(cluster.moduleuse);
				fprintf(fid,'module use %s\n',cluster.moduleuse{i});
			end
			for i=1:length(cluster.moduleload);
				fprintf(fid,'module load %s\n',cluster.moduleload{i});
			end

			fprintf(fid,'export ISSM_DIR="%s/../"\n',cluster.codepath); %FIXME
			fprintf(fid,'source $ISSM_DIR/etc/environment.sh\n');       %FIXME

			fprintf(fid,'\n# Switch to run directory (if not using -l wd):\n');
			fprintf(fid,'cd %s/%s\n\n',cluster.executionpath,dirname);

			if ~md.debug.valgrind
				fprintf(fid,'mpiexec -np %i %s/%s %s %s %s\n',cluster.nprocs(),cluster.codepath,executable,solution,[cluster.executionpath '/' dirname],modelname);
			else
				% Example valgrind usage
				supstring = '';
            % If you have self.valgrindsup or so, otherwise remove
            % for supfile in self.valgrindsup:
            %     supstring += ' --suppressions=' + supfile

				fprintf(fid,'mpiexec -np %i --leak-check=full%s %s/%s %s %s %s\n',cluster.nprocs(),supstring,cluster.codepath,executable,solution,[cluster.executionpath '/' dirname],modelname);
			end

			if ~io_gather, %concatenate the output files:
				fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
			end
			fclose(fid);
		end %}}}
		function UploadQueueJob(cluster,modelname,dirname,filelist) % {{{
			cluster_defaults.UploadQueueJob(cluster,modelname,dirname,filelist);
		end %}}}
		function LaunchQueueJob(cluster,modelname,dirname,filelist,restart,batch) % {{{
			if ~isempty(restart)
				launchcommand=['cd ' cluster.executionpath ' && cd ' dirname ' && hostname && qsub ' modelname '.queue '];
			else
				if strcmpi(cluster.login,oshostname)
					launchcommand=['cd ' cluster.executionpath ' && rm -rf ./' dirname ' && mkdir ' dirname ...
						' && cd ' dirname ' && mv ../' dirname '.tar.gz ./ && tar -zxf ' dirname '.tar.gz  && hostname && qsub ' modelname '.queue '];
				else
					launchcommand=['cd ' cluster.executionpath ' && cd ' dirname ' && hostname && qsub ' modelname '.queue '];
				end
			end
			issmssh(cluster.name,cluster.login,cluster.port,launchcommand);
			%cluster_defaults.LaunchQueueJobSbatch(cluster, modelname, dirname, filelist, restart, batch, 3);
		end %}}}
		function Download(cluster,dirname,filelist) % {{{
			cluster_defaults.Download(cluster,dirname,filelist);
		end %}}}
	end
end
