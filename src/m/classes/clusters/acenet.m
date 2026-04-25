%ACENET class definition
%
%   Usage:
%      cluster=acenet();
%      cluster=acenet('np',3);
%      cluster=acenet('np',3,'login','username');

classdef acenet
	properties (SetAccess=public)
		% {{{
		name='placentia.ace-net.ca'
		login='';
		np=10;
		port=0;
		queue='longq';
		time=10;
		codepath='';
		executionpath='';
		%}}}
	end
	methods
		function cluster=acenet(varargin) % {{{
			%use provided options to change fields
			options=pairoptions(varargin{:});
			%initialize cluster using user settings if provided
			if (exist([cluster.name '_settings'])==2), eval([cluster.name '_settings']); end

			%OK get other fields
			cluster=AssignObjectFields(pairoptions(varargin{:}),cluster);
		end
		%}}} 
		function disp(cluster) % {{{
			% display the object
			disp(sprintf('class ''%s'' object ''%s'' = ',class(cluster),inputname(1)));
			disp(sprintf('    name: %s',cluster.name));
			disp(sprintf('    login: %s',cluster.login));
			disp(sprintf('    np: %i',cluster.np));
			disp(sprintf('    port: %i',cluster.port));
			disp(sprintf('    queue: %s',cluster.queue));
			disp(sprintf('    time: %i',cluster.time));
			disp(sprintf('    codepath: %s',cluster.codepath));
			disp(sprintf('    executionpath: %s',cluster.executionpath));
		end
		%}}}
		function md = checkconsistency(cluster,md,solution,analyses) % {{{

			available_queues={'debug','shortq','longq'};
			queue_requirements_time=[48*1 48*7 48*15];
			queue_requirements_np=[32 128 256];
			QueueRequirements(available_queues,queue_requirements_time,queue_requirements_np,cluster.queue,cluster.np,cluster.time)
		end
		%}}}
		function BuildQueueScript(cluster, md, filename) % {{{

			%Get variables from md
			dirname         = md.private.runtimename;
			modelname       = md.miscellaneous.name;
			solution        = md.private.solution;
			io_gather       = md.settings.io_gather;
			isvalgrind      = md.debug.valgrind;
			isgprof         = md.debug.gprof;
			isdakota        = md.qmu.isdakota;
			isoceancoupling = md.transient.isoceancoupling;

         %checks
			if(isvalgrind) disp('valgrind not supported by cluster, ignoring...'); end
			if(isgprof)    disp('gprof not supported by cluster, ignoring...'); end

			%write queuing script 
			fid=fopen(filename, 'w');
			fprintf(fid,'#!/bin/bash\n');
			fprintf(fid,'#$ -cwd\n');

			fprintf(fid,'#$ -N issm\n');
			tstr = sprintf('#$ -l h_rt=%i:00:00\n',cluster.time);
			fprintf(fid,tstr);

			fprintf(fid,'#$ -l h_vmem=2G\n');

			if strcmp(cluster.executionpath,'/home/klemorza/scratch/issmres.dir')
				% ---- Which acent queue to use ----
				fprintf(fid,'#$ -q short.q@*,medium.q@*\n');
				if cluster.np==4
					% -------- All cpus in the same node --------          
					fprintf(fid,'#$ -pe openmp %i\n',cluster.np);
				else
					fprintf(fid,'#$ -pe ompi %i\n',cluster.np); % To avoid green acenet that does not have InfiniBand
				end

			elseif strcmp(cluster.executionpath,'/net/glacdyn-data/glacdyn/1/klemorza/issm.dir')
				fprintf(fid,'#$ -q tarasov.q\n');
				fprintf(fid,'#$ -l h=cl27*|cl28*|cl29*|cl30*|cl31*|cl320|cl267|cl268|cl269|cl338 \n');
				
				if cluster.np==4
					% -------- All cpus in the same node --------          
					fprintf(fid,'#$ -pe openmp %i\n',cluster.np);
				else
					fprintf(fid,'#$ -pe ompi* %i\n',cluster.np);
				end
			end
						
			fprintf(fid,'#$ -j y\n');
			fprintf(fid,'module purge\n');
			fprintf(fid,'module load intel/12.1.7.367\n');
			fprintf(fid,'module load openmpi/intel/1.2.9\n');
			fprintf(fid,'module load gsl\n');
			fprintf(fid,'export ISSM_DIR="%s/../"\n',cluster.codepath); %FIXME
			fprintf(fid,'source $ISSM_DIR/etc/environment.sh\n');       %FIXME
			fprintf(fid,'\n');
			fprintf(fid,'mpiexec %s/issm.exe %s %s %s 2> %s.errlog >%s.outlog\n',...
					cluster.codepath,solution,[cluster.executionpath '/' dirname],modelname,modelname,modelname);
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
		function LaunchQueueJob(cluster,modelname,dirname,filelist,restart) % {{{

			%Execute Queue job
			if ~isempty(restart)
				launchcommand=['cd ' cluster.executionpath ' && cd ' dirname ' && qsub ' modelname '.queue '];
			else
				launchcommand=['cd ' cluster.executionpath ' && rm -rf ./' dirname ' && mkdir ' dirname ...
					' && cd ' dirname ' && mv ../' dirname '.tar.gz ./ && tar -zxf ' dirname '.tar.gz && qsub ' modelname '.queue '];
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
