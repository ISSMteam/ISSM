%RAIJIN class definition
%
%   Usage:
%      cluster=raijin();
%      cluster=raijin('np',3);
%      cluster=raijin('np',3,'login','username');

classdef raijin
	properties (SetAccess=public)
		% {{{
		name='';
		login='';
		np=0;
		mem=0;
		time=0;
		project='';
		email='';
		port=0;
		queue='';
		codepath='';
		executionpath='';
		modules={};
	end
	%}}}
	methods
		function cluster=raijin(varargin) % {{{

			%initialize cluster using default settings if provided
			if (exist('raijin_settings')==2), raijin_settings; end

			%use provided options to change fields
			cluster=AssignObjectFields(pairoptions(varargin{:}),cluster);
		end
		%}}}
		function disp(cluster) % {{{
			%  display the object
			disp(sprintf('class ''%s'' object ''%s'' = ',class(cluster),inputname(1)));
			disp(sprintf('    name: %s',cluster.name));
			disp(sprintf('    login: %s',cluster.login));
			disp(sprintf('    np: %i',cluster.np));
			disp(sprintf('    mem: %i',cluster.mem));
			disp(sprintf('    time: %i',cluster.time));
			disp(sprintf('    project: %s',cluster.project));
			disp(sprintf('    email: %s',cluster.email));
			disp(sprintf('    port: %i',cluster.port));
			disp(sprintf('    queue: %s',cluster.queue));
			disp(sprintf('    codepath: %s',cluster.codepath));
			disp(sprintf('    executionpath: %s',cluster.executionpath));
			disp(sprintf('    modules: %s',strjoin(cluster.modules,', ')));
		end
		%}}}
		function md = checkconsistency(cluster,md,solution,analyses) % {{{

			available_queues={'express', 'normal'};
			queue_requirements_time=[5*60 40*60];
			queue_requirements_np=[256 1024];

			QueueRequirements(available_queues,queue_requirements_time,queue_requirements_np,cluster.queue,cluster.np,1)

			%Miscellaneous
			if isempty(cluster.login), md = checkmessage(md,'login empty'); end
			if isempty(cluster.codepath), md = checkmessage(md,'codepath empty'); end
			if isempty(cluster.executionpath), md = checkmessage(md,'executionpath empty'); end

		end
		%}}}
		function BuildQueueScript(cluster,dirname,modelname,solution,io_gather,isvalgrind,isgprof,isdakota,isoceancoupling) % {{{

			if(isvalgrind), disp('valgrind not supported by cluster, ignoring...'); end
			if(isgprof),    disp('gprof not supported by cluster, ignoring...'); end

			%write queuing script
			fid=fopen([modelname '.queue'],'w');
			fprintf(fid,'#PBS -S /bin/bash\n');
			fprintf(fid,'#PBS -P %s\n', cluster.project);
			fprintf(fid,'#PBS -q %s\n',cluster.queue);
			fprintf(fid,'#PBS -l ncpus=%i\n',cluster.np);
			fprintf(fid,'#PBS -l walltime=%i\n',cluster.time*60); % walltime is in seconds
			fprintf(fid,'#PBS -l mem=%igb\n', cluster.mem);
			fprintf(fid,'#PBS -M %s\n', cluster.email);
			fprintf(fid,'#PBS -o %s.outlog\n',(modelname));
			fprintf(fid,'#PBS -e %s.errlog\n',(modelname));
			fprintf(fid,'#PBS -l wd\n\n');
			fprintf(fid,'source ~/.bashrc\n');
			fprintf(fid,'source $ISSM_DIR/etc/environment.sh\n');
			for i=1:numel(cluster.modules)
				fprintf(fid,'module load %s\n', cluster.modules{i});
			end
			fprintf(fid,'mpiexec -np %i %s/issm.exe %s %s %s\n',cluster.np,cluster.codepath,solution,[cluster.executionpath '/' dirname],modelname);
			if ~io_gather, %concatenate the output files:
				fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
			end
			fclose(fid);
		end %}}}
		function UploadQueueJob(cluster,modelname,dirname,filelist) % {{{

			%compress the files into one zip.
			compressstring=['tar -zcf ' dirname '.tar.gz '];
			for i=1:numel(filelist),
				compressstring = [compressstring ' ' filelist{i}];
			end
			system(compressstring);

			disp('uploading input file and queuing script');
			issmscpout(cluster.name,cluster.executionpath,cluster.login,cluster.port,{[dirname '.tar.gz']});

		end %}}}
		function LaunchQueueJob(cluster,modelname,dirname,filelist,restart,batch) % {{{
			
			disp('launching solution sequence on remote cluster');
			if ~isempty(restart)
				launchcommand=['cd ' cluster.executionpath ' && cd ' dirname ' && hostname && qsub ' modelname '.queue '];
			else
				launchcommand=['cd ' cluster.executionpath ' && rm -rf ./' dirname ' && mkdir ' dirname ...
					' && cd ' dirname ' && mv ../' dirname '.tar.gz ./ && tar -zxf ' dirname '.tar.gz  && hostname && qsub ' modelname '.queue '];
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
