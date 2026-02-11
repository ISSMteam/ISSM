%PFE cluster class definition
%
%   Usage:
%      cluster=pfe();
%      cluster=pfe('np',3);
%      cluster=pfe('np',3,'login','username');

classdef pfe
	properties (SetAccess=public)
		% {{{
		name           = 'pfe'
		login          = '';
		modules        = {'comp-intel/2018.3.222' '/nasa/intel/impi/2021.3/modulefiles/mpi/2021.3.0' 'scicon/app-tools'};
		numnodes       = 20;
		cpuspernode    = 8;
		port           = 1025;
		queue          = 'long';
		time           = 12*60;
		processor      = 'ivy';
		srcpath        = '';
		extpkgpath     = '';
		codepath       = '';
		executionpath  = '';
		grouplist      = '';
		interactive    = 0;
		bbftp          = 0;
		numstreams     = 8;
		hyperthreading = 0;
	end
	%}}}
	methods
		function cluster=pfe(varargin) % {{{

			%initialize cluster using default settings if provided
			if (exist('pfe_settings')==2), pfe_settings; end

			%use provided options to change fields
			cluster=AssignObjectFields(pairoptions(varargin{:}),cluster);
		end
		%}}}
		function disp(cluster) % {{{
			% display the object
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
			disp(sprintf('    extpkgpath: %s',cluster.extpkgpath));
			disp(sprintf('    codepath: %s',cluster.codepath));
			disp(sprintf('    executionpath: %s',cluster.executionpath));
			disp(sprintf('    grouplist: %s',cluster.grouplist));
			disp(sprintf('    interactive: %i',cluster.interactive));
			disp(sprintf('    bbftp: %s',cluster.bbftp));
			disp(sprintf('    numstreams: %s',cluster.numstreams));
			disp(sprintf('    hyperthreading: %s',cluster.hyperthreading));
		end
		%}}}
		function numprocs=nprocs(cluster) % {{{
			%compute number of processors
			numprocs=cluster.numnodes*cluster.cpuspernode;
		end
		%}}}
		function md = checkconsistency(cluster,md,solution,analyses) % {{{

			available_queues={'long','normal','debug','devel','alphatst@pbspl233'};
			queue_requirements_time=[5*24*60 8*60 2*60 2*60 24*60];
			queue_requirements_np=[2048 2048 150 150 2048];

			QueueRequirements(available_queues,queue_requirements_time,queue_requirements_np,cluster.queue,cluster.nprocs(),cluster.time)

			%now, check cluster.cpuspernode according to processor type
			if strcmpi(cluster.processor,'ivy'),
				if cluster.hyperthreading,
					if ((cluster.cpuspernode>40 ) | (cluster.cpuspernode<1)),
						md = checkmessage(md,'cpuspernode should be between 1 and 40 for ''ivy'' processors in hyperthreading mode');
					end
				else
					if ((cluster.cpuspernode>20 ) | (cluster.cpuspernode<1)),
						md = checkmessage(md,'cpuspernode should be between 1 and 20 for ''ivy'' processors');
					end
				end
			elseif strcmpi(cluster.processor,'bro'),
				if cluster.hyperthreading,
					if ((cluster.cpuspernode>56 ) | (cluster.cpuspernode<1)),
						md = checkmessage(md,'cpuspernode should be between 1 and 56 for ''bro'' processors in hyperthreading mode');
					end
				else
					if ((cluster.cpuspernode>28 ) | (cluster.cpuspernode<1)),
						md = checkmessage(md,'cpuspernode should be between 1 and 28 for ''bro'' processors');
					end
				end
			elseif strcmpi(cluster.processor,'has'),
				if cluster.hyperthreading,
					if ((cluster.cpuspernode>48 ) | (cluster.cpuspernode<1)),
						md = checkmessage(md,'cpuspernode should be between 1 and 48 for ''has'' processors in hyperthreading mode');
					end
				else
					if ((cluster.cpuspernode>24 ) | (cluster.cpuspernode<1)),
						md = checkmessage(md,'cpuspernode should be between 1 and 24 for ''has'' processors');
					end
				end
			
			elseif strcmpi(cluster.processor,'san'),
				if cluster.hyperthreading,
					if ((cluster.cpuspernode>32 ) | (cluster.cpuspernode<1)),
						md = checkmessage(md,'cpuspernode should be between 1 and 32 for ''san'' processors in hyperthreading mode');
					end
				else
					if ((cluster.cpuspernode>16 ) | (cluster.cpuspernode<1)),
						md = checkmessage(md,'cpuspernode should be between 1 and 16 for ''san'' processors');
					end
				end

			elseif strcmpi(cluster.processor,'cas_ait'),
				if cluster.hyperthreading,
					if ((cluster.cpuspernode>80 ) | (cluster.cpuspernode<1)),
						md = checkmessage(md,'cpuspernode should be between 1 and 80 for ''cas_ait'' processors in hyperthreading mode');
					end
				else
					if ((cluster.cpuspernode>40 ) | (cluster.cpuspernode<1)),
						md = checkmessage(md,'cpuspernode should be between 1 and 40 for ''cas_ait'' processors');
					end
				end
			
			elseif strcmpi(cluster.processor,'rom_ait'),
				if cluster.hyperthreading,
					if ((cluster.cpuspernode>256 ) | (cluster.cpuspernode<1)),
						md = checkmessage(md,'cpuspernode should be between 1 and 80 for ''rom_ait'' processors in hyperthreading mode');
					end
				else
					if ((cluster.cpuspernode>128 ) | (cluster.cpuspernode<1)),
						md = checkmessage(md,'cpuspernode should be between 1 and 128 for ''rom_ait'' processors');
					end
				end
			
			else
				md = checkmessage(md,'unknown processor type, should be ''bro'', ''has'', ''ivy'', ''san'', ''cas_ait'', or ''rom_ait''');
			end

			%Miscellaneous
			if isempty(cluster.login), md = checkmessage(md,'login empty'); end
			if isempty(cluster.srcpath), md = checkmessage(md,'srcpath empty'); end
			if isempty(cluster.codepath), md = checkmessage(md,'codepath empty'); end
			if isempty(cluster.executionpath), md = checkmessage(md,'executionpath empty'); end
			if isempty(cluster.grouplist), md = checkmessage(md,'grouplist empty'); end

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
			fprintf(fid,'#PBS -S /bin/bash\n');
%			fprintf(fid,'#PBS -N %s\n',modelname);
			fprintf(fid,'#PBS -l select=%i:ncpus=%i:model=%s\n',cluster.numnodes,cluster.cpuspernode,cluster.processor);
			fprintf(fid,'#PBS -l walltime=%i\n',cluster.time*60); %walltime is in seconds.
			fprintf(fid,'#PBS -q %s\n',cluster.queue);
			fprintf(fid,'#PBS -W group_list=%s\n',cluster.grouplist);
			fprintf(fid,'#PBS -m e\n');
			fprintf(fid,'#PBS -o %s/%s/%s.outlog \n',cluster.executionpath,dirname,modelname);
			fprintf(fid,'#PBS -e %s/%s/%s.errlog \n\n',cluster.executionpath,dirname,modelname);
			fprintf(fid,'. /usr/share/modules/init/bash\n\n');
			for i=1:numel(cluster.modules), fprintf(fid,['module load ' cluster.modules{i} '\n']); end
			fprintf(fid,'export PATH="$PATH:."\n\n');
			fprintf(fid,'export MPI_LAUNCH_TIMEOUT=520\n');
			fprintf(fid,'export MPI_GROUP_MAX=64\n\n');
			fprintf(fid,'export ISSM_DIR="%s"\n',cluster.srcpath); %FIXME
			if cluster.extpkgpath
				fprintf(fid,'export ISSM_EXT_DIR="%s"\n',cluster.extpkgpath); 
			end
			fprintf(fid,'source $ISSM_DIR/etc/environment.sh\n');       %FIXME
			fprintf(fid,'cd %s/%s/\n\n',cluster.executionpath,dirname);
			if ~isvalgrind,
				%fprintf(fid,'/u/scicon/tools/bin/several_tries mpiexec -np %i mbind.x -cs -n%i %s/%s %s %s/%s %s\n',cluster.nprocs(),cluster.cpuspernode,cluster.codepath,executable,solution,cluster.executionpath,dirname,modelname);
				fprintf(fid,'mpiexec -np %i %s/%s %s %s/%s %s\n',cluster.nprocs(),cluster.codepath,executable,solution,cluster.executionpath,dirname,modelname);
			else
				fprintf(fid,'mpiexec -np %i valgrind --leak-check=full %s/%s %s %s %s\n',cluster.nprocs(),cluster.codepath,executable,solution,[cluster.executionpath '/' dirname],modelname);
			end
			if ~io_gather, %concatenate the output files:
				fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
			end
			fclose(fid);

			%in interactive mode, create a run file, and errlog and outlog file
			if cluster.interactive,
				fid=fopen([modelname '.run'],'w');
				if cluster.interactive==10,
						fprintf(fid,'module unload mpi-mvapich2/1.4.1/gcc\n');
						fprintf(fid,'mpiexec -np %i %s/%s %s %s %s\n',cluster.nprocs(),cluster.codepath,executable,solution,[pwd() '/run'],modelname);
				else
					if ~isvalgrind,
						fprintf(fid,'mpiexec -np %i %s/%s %s %s %s\n',cluster.nprocs(),cluster.codepath,executable,solution,[cluster.executionpath '/Interactive' num2str(cluster.interactive)],modelname);
					else
						fprintf(fid,'mpiexec -np %i valgrind --leak-check=full %s/%s %s %s %s\n',cluster.nprocs(),cluster.codepath,executable,solution,[cluster.executionpath '/Interactive' num2str(cluster.interactive)],modelname);
					end
				end
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
		function BuildQueueScriptMultipleModels(cluster,dirname,modelname,solution,dirnames,modelnames,nps) % {{{

			%some checks: 
			if isempty(modelname), error('BuildQueueScriptMultipleModels error message: need a non empty model name!');end

			%what is the executable being called? 
			executable='issm_slc.exe';

			if ispc & ~ismingw, error('BuildQueueScriptMultipleModels not support yet on windows machines');end;

			%write queuing script 
			fid=fopen([modelname '.queue'],'w');

			fprintf(fid,'#PBS -S /bin/bash\n');
			fprintf(fid,'#PBS -N %s\n',modelname);
			fprintf(fid,'#PBS -l select=%i:ncpus=%i:model=%s\n',cluster.numnodes,cluster.cpuspernode,cluster.processor);
			fprintf(fid,'#PBS -l walltime=%i\n',cluster.time*60); %walltime is in seconds.
			fprintf(fid,'#PBS -q %s \n',cluster.queue);
			fprintf(fid,'#PBS -W group_list=%s\n',cluster.grouplist);
			fprintf(fid,'#PBS -m e\n');
			fprintf(fid,'#PBS -o %s.outlog \n',[cluster.executionpath '/' dirname '/' modelname]);
			fprintf(fid,'#PBS -e %s.errlog \n\n',[cluster.executionpath '/' dirname '/' modelname]);
			fprintf(fid,'. /usr/share/modules/init/bash\n\n');
			for i=1:numel(cluster.modules), fprintf(fid,['module load ' cluster.modules{i} '\n']); end
			fprintf(fid,'export PATH="$PATH:."\n\n');
			fprintf(fid,'export MPI_GROUP_MAX=64\n\n');
			fprintf(fid,'export ISSM_DIR="%s/../"\n',cluster.codepath); %FIXME
			fprintf(fid,'source $ISSM_DIR/etc/environment.sh\n');       %FIXME
			fprintf(fid,'cd %s/%s/\n\n',cluster.executionpath,dirname);

			%number of cpus: 
			mpistring=sprintf('mpiexec -np %i ',cluster.numnodes*cluster.cpuspernode);

			%executable: 
			mpistring=[mpistring sprintf('%s/%s ',cluster.codepath,executable)];

			%solution name: 
			mpistring=[mpistring sprintf('%s ',solution)];

			%execution directory and model name: 
			mpistring=[mpistring sprintf('%s/%s %s',cluster.executionpath,dirname,modelname)];

			%inform main executable of how many icecaps, glaciers and earth models are being run: 
			mpistring=[mpistring sprintf(' %i ',length(dirnames))];

			%icecaps, glaciers and earth location, names and number of processors associated:
			for i=1:length(dirnames),
				mpistring=[mpistring sprintf(' %s/%s %s %i ',cluster.executionpath,dirnames{i},modelnames{i},nps{i})];
			end

			%write this long string to disk: 
			fprintf(fid,mpistring);
			fclose(fid);

			if cluster.interactive,
				fid=fopen([modelname '.run'],'w');

				%number of cpus: 
				mpistring=sprintf('mpiexec -np %i ',cluster.numnodes*cluster.cpuspernode);

				%executable: 
				mpistring=[mpistring sprintf('%s/%s ',cluster.codepath,executable)];

				%solution name: 
				mpistring=[mpistring sprintf('%s ',solution)];

				%execution directory and model name: 
				mpistring=[mpistring sprintf('%s/%s %s',cluster.executionpath,dirname,modelname)];

				%inform main executable of how many icecaps, glaciers and earth models are being run: 
				mpistring=[mpistring sprintf(' %i ',length(dirnames))];

				%icecaps, glaciers and earth location, names and number of processors associated:
				for i=1:length(dirnames),
					mpistring=[mpistring sprintf(' %s/Interactive%i %s %i ',cluster.executionpath,cluster.interactive,modelnames{i},nps{i})];
				end

				%write this long string to disk: 
				fprintf(fid,mpistring);
				fclose(fid);

				fid=fopen([modelname '.errlog'],'w');
				fclose(fid);
				fid=fopen([modelname '.outlog'],'w');
				fclose(fid);
			end
		end
		%}}}
		function BuildKrigingQueueScript(cluster,modelname,solution,io_gather,isvalgrind,isgprof) % {{{

			if(isgprof),    disp('gprof not supported by cluster, ignoring...'); end

			%write queuing script 
			fid=fopen([modelname '.queue'],'w');
			fprintf(fid,'#PBS -S /bin/bash\n');
			%			fprintf(fid,'#PBS -N %s\n',modelname);
			fprintf(fid,'#PBS -l select=%i:ncpus=%i:model=%s\n',cluster.numnodes,cluster.cpuspernode,cluster.processor);
			fprintf(fid,'#PBS -l walltime=%i\n',cluster.time*60); %walltime is in seconds.
			fprintf(fid,'#PBS -q %s \n',cluster.queue);
			fprintf(fid,'#PBS -W group_list=%s\n',cluster.grouplist);
			fprintf(fid,'#PBS -m e\n');
			fprintf(fid,'#PBS -o %s.outlog \n',modelname);
			fprintf(fid,'#PBS -e %s.errlog \n\n',modelname);
			fprintf(fid,'. /usr/share/modules/init/bash\n\n');
			for i=1:numel(cluster.modules), fprintf(fid,['module load ' cluster.modules{i} '\n']); end
			fprintf(fid,'export PATH="$PATH:."\n');
			fprintf(fid,'export ISSM_DIR="%s/../"\n',cluster.codepath); %FIXME
			fprintf(fid,'source $ISSM_DIR/etc/environment.sh\n');       %FIXME
			fprintf(fid,'export MPI_GROUP_MAX=64\n\n');
			fprintf(fid,'cd %s/%s/\n\n',cluster.executionpath,modelname);
			fprintf(fid,'mpiexec -np %i %s/kriging.exe %s %s\n',cluster.nprocs(),cluster.codepath,[cluster.executionpath '/' modelname],modelname); %FIXME
			if ~io_gather, %concatenate the output files:
				fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
			end
			fclose(fid);

			%in interactive mode, create a run file, and errlog and outlog file
			if cluster.interactive,
				fid=fopen([modelname '.run'],'w');
				if ~isvalgrind,
					fprintf(fid,'mpiexec -np %i %s/kriging.exe %s %s\n',cluster.nprocs(),cluster.codepath,[cluster.executionpath '/' modelname],modelname);
				else
					fprintf(fid,'mpiexec -np %i valgrind --leak-check=full %s/kriging.exe %s %s\n',cluster.nprocs(),cluster.codepath,[cluster.executionpath '/' modelname],modelname);
				end
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
		function BuildOceanQueueScript(np,cluster,modelname) % {{{

			%write queuing script 
			fid=fopen([modelname '.queue'],'w');
			fprintf(fid,'#PBS -S /bin/bash\n');
			fprintf(fid,'#PBS -l select=1:ncpus=%i:model=%s\n',np,cluster.processor);
			fprintf(fid,'#PBS -l walltime=%i\n',cluster.time*60); %walltime is in seconds.
			fprintf(fid,'#PBS -q %s \n',cluster.queue);
			fprintf(fid,'#PBS -W group_list=%s\n',cluster.grouplist);
			fprintf(fid,'#PBS -m e\n');
			fprintf(fid,'#PBS -o %s.outlog \n',modelname);
			fprintf(fid,'#PBS -e %s.errlog \n\n',modelname);
			fprintf(fid,'. /usr/share/modules/init/bash\n\n');
			%for i=1:numel(cluster.modules), fprintf(fid,['module load ' cluster.modules{i} '\n']); end %FIXME: should use this!
			fprintf(fid,'module load comp-intel/2016.2.181\n');
			fprintf(fid,'module load netcdf/4.4.1.1_mpt\n');
			fprintf(fid,'module load mpi-sgi/mpt.2.15r20\n');
			fprintf(fid,'export PATH="$PATH:."\n');
			fprintf(fid,'export MPI_GROUP_MAX=64\n\n');
			fprintf(fid,['cd ' pwd() ' \n\n']);
			fprintf(fid,'mpiexec -np %i ./mitgcmuv\n',np); 
		%	if ~io_gather, %concatenate the output files:
		%		fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
		%	end
			fclose(fid);

			%in interactive mode, create a run file, and errlog and outlog file
			if cluster.interactive,
				fid=fopen([modelname '.run'],'w');
				fprintf(fid,'module load mpi-sgi/mpt.2.15r20\n');
				fprintf(fid,['mpiexec -np %i ./mitgcmuv \n'],np);
				fprintf(fid,['touch ' modelname '.lock %s\n']);
				fclose(fid);
				fid=fopen([modelname '.errlog'],'w');
				fclose(fid);
				fid=fopen([modelname '.outlog'],'w');
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

			disp('uploading input file and queueing script');
			if cluster.interactive==10,
				directory=[pwd() '/run/'];
			elseif cluster.interactive,
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

			%launch command, to be executed via ssh
			if cluster.interactive,
				if ~isempty(restart)
					launchcommand=['cd ' cluster.executionpath '/Interactive' num2str(cluster.interactive)];
				else
					if cluster.interactive==10,
						launchcommand=['cd ' pwd() '/run && tar -zxf ' dirname '.tar.gz'];
					else
						launchcommand=['cd ' cluster.executionpath '/Interactive' num2str(cluster.interactive) ' && tar -zxf ' dirname '.tar.gz'];
					end
				end
			else
				if ~isempty(restart)
					launchcommand=['cd ' cluster.executionpath ' && cd ' dirname ' && /PBS/bin/qsub ' modelname '.queue '];
				else
					launchcommand=['cd ' cluster.executionpath ' && rm -rf ./' dirname ' && mkdir ' dirname ...
						' && cd ' dirname ' && mv ../' dirname '.tar.gz ./ && tar -zxf ' dirname '.tar.gz && /PBS/bin/qsub ' modelname '.queue '];
				end
			end

			%Execute Queue job
			issmssh(cluster.name,cluster.login,cluster.port,launchcommand);
		end
		%}}}
		function Download(cluster,dirname,filelist) % {{{

			%copy files from cluster to current directory
			if cluster.interactive==10,
				directory=[pwd() '/run/'];
			elseif ~cluster.interactive,
				directory=[cluster.executionpath '/' dirname '/'];
			else
				directory=[cluster.executionpath '/Interactive' num2str(cluster.interactive) '/'];
			end

			if cluster.bbftp,
				issmbbftpin(cluster.name, cluster.login, cluster.port, cluster.numstreams, directory, filelist);
			else
				issmscpin(cluster.name,cluster.login,cluster.port,directory,filelist);
			end

		end %}}}
	end
end
