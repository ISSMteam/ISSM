%AWS_ISSM_SOLUTION_SERVER cluster class definition
%
%   Usage:
%      cluster=aws_issm_solution_server();
%      cluster=aws_issm_solution_server('np',3);
%      cluster=aws_issm_solution_server('np',3,'login','username');
%
%   TODO:
%   - Modify src/m/os/issmscp* with idfile parameter

classdef aws_issm_solution_server
	properties (SetAccess=public)
		% {{{
		name           = '54.67.123.214';
		login          = '';
		idfile         = '';
		modules        = {};
		numnodes       = 1;
		cpuspernode    = 8;
		queue          = '';
		time           = 12*60*60;
		srcpath        = '/usr/local/repos/issm/trunk-jpl-head';
		extpkgpath     = '/usr/local/issm-ext';
		codepath       = '/usr/local/repos/issm/trunk-jpl-head/bin';
		executionpath  = '';
		interactive    = 0;
		numstreams     = 1;
		hyperthreading = 0;
		email          = '';
	end
	%}}}
	methods
		function cluster=aws_issm_solution_server(varargin) % {{{

			%initialize cluster using default settings if provided
			if (exist('aws_issm_solution_server_settings')==2), aws_issm_solution_server_settings; end

			%use provided options to change fields
			cluster=AssignObjectFields(pairoptions(varargin{:}),cluster);
		end
		%}}}
		function disp(cluster) % {{{
			% display the object
			disp(sprintf('class ''%s'' object ''%s'' = ',class(cluster),inputname(1)));
			disp(sprintf('    name: %s',cluster.name));
			disp(sprintf('    login: %s',cluster.login));
			disp(sprintf('    idfile: %s',cluster.idfile));
			disp(sprintf('    modules: %s',strjoin(cluster.modules,', ')));
			disp(sprintf('    numnodes: %i',cluster.numnodes));
			disp(sprintf('    cpuspernode: %i',cluster.cpuspernode));
			disp(sprintf('    np: %i',cluster.nprocs()));
			disp(sprintf('    time: %i',cluster.time));
			disp(sprintf('    processor: %i',cluster.processor));
			disp(sprintf('    srcpath: %s',cluster.srcpath));
			disp(sprintf('    extpkgpath: %s',cluster.extpkgpath));
			disp(sprintf('    codepath: %s',cluster.codepath));
			disp(sprintf('    executionpath: %s',cluster.executionpath));
			disp(sprintf('    grouplist: %s',cluster.grouplist));
			disp(sprintf('    interactive: %i',cluster.interactive));
			disp(sprintf('    numstreams: %s',cluster.numstreams));
			disp(sprintf('    hyperthreading: %s',cluster.hyperthreading));
			disp(sprintf('    email: %s',cluster.email));
		end
		%}}}
		function numprocs=nprocs(cluster) % {{{
			%compute number of processors
			numprocs=cluster.numnodes*cluster.cpuspernode;
		end
		%}}}
		function md = checkconsistency(cluster,md,solution,analyses) % {{{

			if ((cluster.numnodes>1 ) | (cluster.numnodes<1)),
				md = checkmessage(md,'only 1 node is currently available');
			end

			if ((cluster.cpuspernode>8 ) | (cluster.cpuspernode<1)),
				md = checkmessage(md,'cpuspernode should be between 1 and 8');
			end

			%Miscellaneous
			if isempty(cluster.login), md = checkmessage(md,'login empty'); end
			if isempty(cluster.srcpath), md = checkmessage(md,'srcpath empty'); end
			if isempty(cluster.extpkgpath), md = checkmessage(md,'extpkgpath empty'); end
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
			fprintf(fid,'export PATH="${PATH}:."\n');
			fprintf(fid,'export MPI_LAUNCH_TIMEOUT=520\n');
			fprintf(fid,'export MPI_GROUP_MAX=64\n');
			fprintf(fid,'export ISSM_DIR="%s"\n',cluster.srcpath);
			if cluster.extpkgpath
				fprintf(fid,'export ISSM_EXT_DIR="%s"\n',cluster.extpkgpath);
			end
			fprintf(fid,'source $ISSM_DIR/etc/environment.sh\n\n');
			if cluster.interactive
				if IssmConfig('_HAVE_MPI_'),
					fprintf(fid,'mpiexec -np %i %s/%s %s %s %s\n',cluster.nprocs(),cluster.codepath,executable,solution,[cluster.executionpath '/' dirname],modelname);
				else
					fprintf(fid,'%s/%s %s %s %s\n',cluster.codepath,executable,solution,[cluster.executionpath '/' dirname],modelname);
				end
			else
				if IssmConfig('_HAVE_MPI_'),
					fprintf(fid,'mpiexec -np %i %s/%s %s %s %s 2> %s.errlog > %s.outlog &\n',cluster.nprocs(),cluster.codepath,executable,solution,[cluster.executionpath '/' dirname],modelname,modelname,modelname);
				else
					fprintf(fid,'%s/%s %s %s %s 2> %s.errlog > %s.outlog &\n',cluster.codepath,executable,solution,[cluster.executionpath '/' dirname],modelname,modelname,modelname);
				end
			end
			if ~io_gather, %concatenate the output files:
				fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
			end
			fclose(fid);

			%in interactive mode, create a run file, and errlog and outlog file
			if cluster.interactive,
				fid=fopen([modelname '.errlog'],'w'); fclose(fid);
				fid=fopen([modelname '.outlog'],'w'); fclose(fid);
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

			%upload input files
			if cluster.interactive==10,
				directory=[pwd() '/run/'];
			elseif cluster.interactive,
				directory=[cluster.executionpath '/Interactive' num2str(cluster.interactive)];
			else 
				directory=cluster.executionpath;
			end

			%NOTE: Replacement for issmscpout(cluster.name,directory,cluster.login,cluster.port,{[dirname '.tar.gz']});
			uploadstring=['scp -i ' cluster.idfile ' ' dirname '.tar.gz ' cluster.login '@' cluster.name ':' directory];
			[status,result]=system(uploadstring);
			if status, 
				error(['cluster.UploadQueueJob error message: ' status]);
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
					launchcommand=['cd ' cluster.executionpath ' && cd ' dirname ' && ./' modelname '.queue'];
				else
					launchcommand=['cd ' cluster.executionpath ' && rm -rf ./' dirname ' && mkdir ' dirname ...
						' && cd ' dirname ' && mv ../' dirname '.tar.gz . && tar -zxf ' dirname '.tar.gz && chmod +x ./' modelname '.queue && ./' modelname '.queue'];
				end
			end

			%Execute Queue job
			%NOTE: Replacement for issmssh(cluster.name,cluster.login,cluster.port,launchcommand);
			launchstring=['ssh -l ' cluster.login ' -i ' cluster.idfile ' ' cluster.name ' "' launchcommand '"'];
			[status,result]=system(launchstring);
			if status,
				error(['cluster.LaunchQueueJob error message: ' status]);
			end
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

			%NOTE: Replacement for issmscpin(cluster.name,cluster.login,cluster.port,directory,filelist);
			if numel(filelist)==1,
				filestring=filelist{1};
			else
				filestring='\{';
				for i=1:numel(filelist)-1,
					filestring=[filestring filelist{i} ','];
				end
				filestring=[filestring filelist{end} '\}'];
			end

			downloadstring=['scp -i ' cluster.idfile ' ' cluster.login '@' cluster.name ':' directory '/' filestring ' ./'];
			[status,result]=system(downloadstring);
			if status, 
				error(['cluster.Download error message: ' status]);
			end

			%check scp worked
			for i=1:numel(filelist),
				if ~exist(['./' filelist{i}]),
					warning(['cluster.Download error message: could not scp ' filelist{i}]);
				end
			end

		end %}}}
	end
end
