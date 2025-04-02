%DISCOVER cluster class definition
%
%   Usage:
%      cluster=discover();
%      cluster=discover('np',3);
%      cluster=discover('np',3,'login','username');

classdef discover
	properties (SetAccess=public)
		% {{{
		name='login.nccs.nasa.gov';
		login='';
		modules={'comp/intel/20.0.0.166' 'mpi/sgi-mpt/2.17' 'cmake/3.17.0'};
		numnodes=20;
		cpuspernode=8;
		port=0;
		queue='general';
		time=12*60*60;
		processor='west';
		codepath='';
		executionpath='';
		grouplist='';
		interactive=0;
		bbftp=0;
		numstreams=8;
		hyperthreading=0;
		email='';
	end
	%}}}
	methods
		function cluster=discover(varargin) % {{{

			%initialize cluster using default settings if provided
			if (exist('discover_settings')==2), discover_settings; end

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
			disp(sprintf('    srcpath: %s',cluster.srcpath));
			disp(sprintf('    codepath: %s',cluster.codepath));
			disp(sprintf('    executionpath: %s',cluster.executionpath));
			disp(sprintf('    grouplist: %s',cluster.grouplist));
			disp(sprintf('    interactive: %i',cluster.interactive));
			disp(sprintf('    bbftp: %s',cluster.bbftp));
			disp(sprintf('    numstreams: %s',cluster.numstreams));
			disp(sprintf('    hyperthreading: %i',cluster.hyperthreading));
			disp(sprintf('    email: %s',cluster.email));
		end
		%}}}
		function numprocs=nprocs(cluster) % {{{
			%compute number of processors
			numprocs=cluster.numnodes*cluster.cpuspernode;
		end
		%}}}
		function md = checkconsistency(cluster,md,solution,analyses) % {{{

			available_queues={'long','allnccs','debug'};
			queue_requirements_time=[24*60*60 12*60*60 1*60*60];
			queue_requirements_np=[560 6000 532];

			QueueRequirements(available_queues,queue_requirements_time,queue_requirements_np,cluster.queue,cluster.nprocs(),cluster.time)

			%now, check cluster.cpuspernode according to processor type
			if ( strcmpi(cluster.processor,'sand')),
				if ((cluster.cpuspernode>16 ) | (cluster.cpuspernode<1)),
					md = checkmessage(md,'cpuspernode should be between 1 and 16 for ''sand'' processors');
				end
			elseif strcmpi(cluster.processor,'hasw'),
				if ((cluster.cpuspernode>28 ) | (cluster.cpuspernode<1)),
					md = checkmessage(md,'cpuspernode should be between 1 and 28 for ''hasw'' processors');
				end
			else
				md = checkmessage(md,'unknown processor type, should be ''sand'' or ''hasw'' ');
			end

			%Miscellaneous
			if isempty(cluster.login), md = checkmessage(md,'login empty'); end
			if isempty(cluster.codepath), md = checkmessage(md,'codepath empty'); end
			if isempty(cluster.executionpath), md = checkmessage(md,'executionpath empty'); end

		end
		%}}}
		function BuildQueueScript(cluster,dirname,modelname,solution,io_gather,isvalgrind,isgprof,isdakota,isoceancoupling) % {{{

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
			fprintf(fid,'#SBATCH --qos=%s \n',cluster.queue);
			fprintf(fid,'#SBATCH -o %s.outlog \n',modelname);
			fprintf(fid,'#SBATCH -e %s.errlog \n',modelname);
			fprintf(fid,'#SBATCH -n %i \n',cluster.nprocs());
			fprintf(fid,'#SBATCH -N %i \n',cluster.numnodes);
			fprintf(fid,'#SBATCH -t %02i:%02i:00 \n',floor(cluster.time/3600),floor(mod(cluster.time,3600)/60));
			fprintf(fid,'#SBATCH -A %s \n\n',cluster.grouplist);
			if length(find(cluster.email=='@'))>0
				fprintf(fid,'#SBATCH --mail-user=%s \n',cluster.email);
				fprintf(fid,'#SBATCH --mail-type=end \n\n');
			end
			fprintf(fid,'. /usr/share/modules/init/bash\n\n');
			for i=1:numel(cluster.modules),
				fprintf(fid,['module load ' cluster.modules{i} '\n']);
			end
			fprintf(fid,'export MPI_GROUP_MAX=64\n\n');
			fprintf(fid,'export MPI_UNBUFFERED_STDIO=true\n\n');
			fprintf(fid,'export PATH="$PATH:."\n\n');
			fprintf(fid,'export ISSM_DIR="%s/../"\n',cluster.codepath); %FIXME
			fprintf(fid,'source $ISSM_DIR/etc/environment.sh\n');       %FIXME
			fprintf(fid,'cd %s/%s\n\n',cluster.executionpath,dirname);

			fprintf(fid,'mpiexec -np %i %s/%s %s %s/%s %s\n',cluster.nprocs(),cluster.codepath,executable,solution,cluster.executionpath,dirname,modelname);
			if ~io_gather, %concatenate the output files:
				fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
			end
			fclose(fid);

			%in interactive mode, create a run file, and errlog and outlog file
			if cluster.interactive,
				fid=fopen([modelname '.run'],'w');
				if ~isvalgrind,
					fprintf(fid,'mpiexec -np %i %s/%s %s %s/%s %s\n',cluster.nprocs(),cluster.codepath,executable,solution,cluster.executionpath,dirname,modelname);
				else
					fprintf(fid,'mpiexec -np %i valgrind --leak-check=full %s/%s %s %s/%s %s\n',cluster.nprocs(),cluster.codepath,executable,solution,cluster.executionpath,dirname,modelname);
				end
				if ~io_gather, %concatenate the output files:
					fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
				end
				fclose(fid);
				fid=fopen([modelname '.errlog'],'w'); % TODO: Change this to system call (touch <file>)?
				fclose(fid);
				fid=fopen([modelname '.outlog'],'w'); % TODO: Change this to system call (touch <file>)?
				fclose(fid);
			end
		end %}}}
		function UploadQueueJob(cluster,modelname,dirname,filelist) % {{{

			%compress the files into one zip.
			compressstring=['tar -zcf ' dirname '.tar.gz'];
			for i=1:numel(filelist),
				compressstring = [compressstring ' ' filelist{i}];
			end
			if cluster.interactive,
				compressstring = [compressstring ' ' modelname '.run '  modelname '.errlog ' modelname '.outlog '];
			end
			system(compressstring);

			disp('uploading input file and queuing script');
			if cluster.interactive,
				directory=[cluster.executionpath '/Interactive' num2str(cluster.interactive)];
			else 
				directory=cluster.executionpath;
			end

			if cluster.bbftp,
				issmbbftpout(cluster.name,directory,cluster.login,cluster.port,cluster.numstreams,{[dirname '.tar.gz']});
			else
				issmscpout(cluster.name,directory,cluster.login,cluster.port,{[dirname '.tar.gz']});
			end

		end
		%}}}
		function LaunchQueueJob(cluster,modelname,dirname,filelist,restart,batch) % {{{

			if cluster.interactive,
				if ~isempty(restart)
					launchcommand=['cd ' cluster.executionpath '/Interactive' num2str(cluster.interactive)];
				else
					launchcommand=['cd ' cluster.executionpath '/Interactive' num2str(cluster.interactive) ' && tar -zxf ' dirname '.tar.gz'];
				end
			else
				if ~isempty(restart)
					launchcommand=['cd ' cluster.executionpath ' && cd ' dirname ' && sbatch ' modelname '.queue '];
				else
					launchcommand=['cd ' cluster.executionpath ' && rm -rf ./' dirname ' && mkdir ' dirname ...
						' && cd ' dirname ' && mv ../' dirname '.tar.gz ./ && tar -zxf ' dirname '.tar.gz  && sbatch ' modelname '.queue '];
				end
			end

			disp('launching solution sequence on remote cluster');
			issmssh(cluster.name,cluster.login,cluster.port,launchcommand);
		end
		%}}}
		function Download(cluster,dirname,filelist) % {{{

			%copy files from cluster to current directory
			if cluster.interactive,
				directory=[cluster.executionpath '/Interactive' num2str(cluster.interactive) '/'];
			else
				directory=[cluster.executionpath '/' dirname '/'];
			end

			if cluster.bbftp,
				issmbbftpin(cluster.name, cluster.login, cluster.port, cluster.numstreams, directory, filelist);
			else
				issmscpin(cluster.name,cluster.login,cluster.port,directory,filelist);
			end

		end %}}}
	end
end
