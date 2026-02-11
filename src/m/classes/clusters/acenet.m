%ACENET class definition
%
%   Usage:
%      cluster=acenet();
%      cluster=acenet('np',3);
%      cluster=acenet('np',3,'login','username');

classdef acenet
	properties (SetAccess=public)
		% {{{
		%name='glacdyn.ace-net.ca'
		name='placentia.ace-net.ca'
		%name='brasdor.ace-net.ca'
		login='klemorza';
		np=10;
		port=0;
		queue='longq';
		time=10;
		% codepath='/usr/local/issm-r11321/bin'; % this one is for issm on acenet global
		codepath='/home/klemorza/issm/trunk-jpl/bin'; % this one is for issm on my acenet directory
		%executionpath='/home/klemorza/issm/trunk-jpl/execution';
		%executionpath='/home/klemorza/scratch/issmres.dir';
		executionpath='/net/glacdyn-data/glacdyn/1/klemorza/issm.dir';
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
		function BuildQueueScript(cluster,dirname,modelname,solution,io_gather,isvalgrind,isgprof,isdakota,isoceancoupling) % {{{

			if(isvalgrind), disp('valgrind not supported by cluster, ignoring...'); end
			if(isgprof),    disp('gprof not supported by cluster, ignoring...'); end

			%write queuing script 
			fid=fopen([modelname '.queue'],'w');
			fprintf(fid,'#!/bin/bash\n');
			fprintf(fid,'#$ -cwd\n');

			fprintf(fid,'#$ -N issm\n');
			% fprintf(fid,'#$ -l h_rt=00:15:00\n');
			% fprintf(fid,'#$ -l h_rt=5:00:0\n');
			% fprintf(fid,'#$ -l h_rt=25:00:0\n');
			% fprintf(fid,'#$ -l h_rt=47:59:00\n');
			% fprintf(fid,'#$ -l h_rt=72:00:0\n');
			% fprintf(fid,'#$ -l h_rt=96:00:0\n');
			% fprintf(fid,'#$ -l h_rt=336:00:0\n');
			tstr = sprintf('#$ -l h_rt=%i:00:00\n',cluster.time);
			fprintf(fid,tstr);

			fprintf(fid,'#$ -l h_vmem=2G\n');

			if strcmp(cluster.executionpath,'/home/klemorza/scratch/issmres.dir')
				% ---- Which acent queue to use ----
				fprintf(fid,'#$ -q short.q@*,medium.q@*\n');
				%fprintf(fid,'#$ -q medium.q@*,long.q@*\n');
				%fprintf(fid,'#$ -q medium.q@*\n');
				%fprintf(fid,'#$ -q short.q@*\n');
				% Acenet nodes with 16cpus and more than 60G mem
				% fprintf(fid,'#$ -l h=cl001|cl002|cl003|cl004|cl005|cl006|cl007|cl008|cl009|cl010|cl011|cl012|cl021|cl022|cl023|cl024 \n');
				% ---- cpus on different nodes ----
				if cluster.np==4
					% -------- All cpus in the same node --------          
					fprintf(fid,'#$ -pe openmp %i\n',cluster.np);
				else
					fprintf(fid,'#$ -pe ompi %i\n',cluster.np); % To avoid green acenet that does not have InfiniBand
				end

			elseif strcmp(cluster.executionpath,'/net/glacdyn-data/glacdyn/1/klemorza/issm.dir')
				% ---- Which node for Lev's queue are selected ----
				fprintf(fid,'#$ -q tarasov.q\n');
				fprintf(fid,'#$ -l h=cl27*|cl28*|cl29*|cl30*|cl31*|cl320|cl267|cl268|cl269|cl338 \n');
				%fprintf(fid,'#$ -l h=cl27*|cl28*|cl29*|cl30*|cl31*|cl320|cl267|cl268|cl269 \n');
				%fprintf(fid,'#$ -l h=cl0* \n');
				% fprintf(fid,'#$ -l h=cl338 \n');
				
				if cluster.np==4
					% -------- All cpus in the same node --------          
					fprintf(fid,'#$ -pe openmp %i\n',cluster.np);
				else
					fprintf(fid,'#$ -pe ompi* %i\n',cluster.np);
					%fprintf(fid,'#$ -pe 4per %i\n',cluster.np);
					%fprintf(fid,'#$ -pe 8per %i\n',cluster.np);
				end
			end
						
			% ---- misc ----
			fprintf(fid,'#$ -j y\n');
			
			fprintf(fid,'module purge\n');
			%fprintf(fid,'module load gcc openmpi/gcc\n');
			%fprintf(fid,'module unload openmpi\n');
			fprintf(fid,'module load intel/12.1.7.367\n');
			fprintf(fid,'module load openmpi/intel/1.2.9\n');

			fprintf(fid,'module load gsl\n');
			%fprintf(fid,'module load issm\n');
			fprintf(fid,'export ISSM_DIR="%s/../"\n',cluster.codepath); %FIXME
			fprintf(fid,'source $ISSM_DIR/etc/environment.sh\n');       %FIXME
			fprintf(fid,'\n');
			fprintf(fid,'mpiexec %s/issm.exe %s %s %s 2> %s.errlog >%s.outlog\n',...
					cluster.codepath,solution,[cluster.executionpath '/' dirname],modelname,modelname,modelname);
			%fprintf(fid,'echo $HOSTNAME >>%s.outlog',modelname);
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
