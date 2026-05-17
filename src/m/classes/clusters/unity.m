%UNITY (Massachusetts Green High Performance Computing Center) cluster class definition
%
%   Usage:
%      cluster=unity();
%      cluster=unity('np',3);
%      cluster=unity('np',3,'login','username');

classdef unity
	properties (SetAccess=public)
		% {{{
		name          = 'unity'
		login         = '';
		numnodes      = 1;
		cpuspernode   = 16;
		codepath      = '';
		executionpath = '';
		time          = 10; %in hours
		memory        = 32;  %in Gb
		email         = '';
	end
	%}}}
	methods
		function cluster=unity(varargin) % {{{

			%initialize cluster using default settings if provided
			if (exist('unity_settings')==2), unity_settings; end

			%use provided options to change fields
			cluster=AssignObjectFields(pairoptions(varargin{:}),cluster);
		end
		%}}}
		function disp(cluster) % {{{
			%  display the object
			disp(sprintf('class ''%s'' object ''%s'' = ',class(cluster),inputname(1)));
			disp(sprintf('    name:  %s',cluster.name));
			disp(sprintf('    login: %s',cluster.login));
			disp(sprintf('    numnodes:    %i',cluster.numnodes));
			disp(sprintf('    cpuspernode: %i',cluster.cpuspernode));
			disp(sprintf('    time: %i hours',cluster.time));
			disp(sprintf('    memory: %i Gb',cluster.memory));
			disp(sprintf('    email: %s (receive notifications if END,FAIL)',cluster.email));
			disp(sprintf('    codepath:      %s',cluster.codepath));
			disp(sprintf('    executionpath: %s',cluster.executionpath));
		end
		%}}}
		function numprocs=nprocs(cluster) % {{{
			%compute number of processors
			numprocs=cluster.numnodes*cluster.cpuspernode;
		end
		%}}}
		function md = checkconsistency(cluster,md,solution,analyses) % {{{
			%Miscellaneous
			if isempty(cluster.login), md = checkmessage(md,'login empty'); end
			if isempty(cluster.codepath), md = checkmessage(md,'codepath empty'); end
			if isempty(cluster.executionpath), md = checkmessage(md,'executionpath empty'); end
		end
		%}}}
		function BuildQueueScript(cluster, md, filename) % {{{

         %Get variables from md
         dirname         = md.private.runtimename;
         modelname       = md.miscellaneous.name;
         solution        = md.private.solution;
         io_gather       = md.settings.io_gather;
         isdakota        = md.qmu.isdakota;
         isoceancoupling = md.transient.isoceancoupling;

			%write queuing script
			fid=fopen(filename, 'w');
			fprintf(fid,'#!/bin/bash -l\n');
			fprintf(fid,'#SBATCH --job-name=%s\n',modelname);
			fprintf(fid,'#SBATCH -p cpu  # Partition\n');
			fprintf(fid,'#SBATCH -o %s.outlog \n',modelname);
			fprintf(fid,'#SBATCH -e %s.errlog \n',modelname);
			fprintf(fid,'#SBATCH --nodes=%i\n',cluster.numnodes);
			fprintf(fid,'#SBATCH --ntasks-per-node=%i\n',cluster.cpuspernode);
			fprintf(fid,'#SBATCH --time=%s\n',eraseBetween(datestr(cluster.time/24,'dd-HH:MM:SS'),1,1)); %walltime is in d-HH:MM:SS format. cluster.time is in hour
			fprintf(fid,'#SBATCH --mem=%iG\n',cluster.memory);
			if ~isempty(cluster.email)
				fprintf(fid,'#SBATCH --mail-type=END,FAIL\n');
				fprintf(fid,'#SBATCH --mail-user=%s\n', cluster.email);
			end
			fprintf(fid,'\n');
			fprintf(fid,'module load intel-oneapi-compilers/2024.1.0 intel-oneapi-mpi/2021.12.1 petsc/3.22.1\n');
			fprintf(fid,'mpiexec -n %i %s/issm.exe %s %s %s\n',cluster.nprocs(), cluster.codepath,solution,[cluster.executionpath '/' dirname],modelname);
			if ~io_gather, %concatenate the output files:
				fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
			end
			fclose(fid);
		end %}}}
		function UploadQueueJob(cluster,modelname,dirname,filelist) % {{{

			%compress the files into one zip.
			%filelist contains full paths; tar with -C so only basenames are stored in the archive
			root=[issmdir() '/execution/' dirname];
			compressstring=['tar -C ' root ' -zcf ' dirname '.tar.gz'];
			for i=1:numel(filelist)
				if ~exist(filelist{i},'file')
					error(['File ' filelist{i} ' not found']);
				end
				[~,fname,fext]=fileparts(filelist{i});
				compressstring=[compressstring ' ' fname fext];
			end
			system(compressstring);

			%upload input files
			issmscpout(cluster.name,cluster.executionpath,cluster.login,0,{[dirname '.tar.gz']});

		end %}}}
		function LaunchQueueJob(cluster,modelname,dirname,filelist,restart,batch) % {{{

			%Execute Queue job
			if ~isempty(restart)
				launchcommand=['cd ' cluster.executionpath ' && cd ' dirname ' && hostname && sbatch ' modelname '.queue '];
			else
				launchcommand=['cd ' cluster.executionpath ' && rm -rf ./' dirname ' && mkdir ' dirname ...
					' && cd ' dirname ' && mv ../' dirname '.tar.gz ./ && tar -zxf ' dirname '.tar.gz && sbatch ' modelname '.queue '];
			end
			issmssh(cluster.name,cluster.login,0,launchcommand);
		end %}}}
		function Download(cluster,dirname,filelist) % {{{

			%copy files from cluster to current directory
			directory=[cluster.executionpath '/' dirname '/'];
			issmscpin(cluster.name,cluster.login,0,directory,filelist, 2); %use {} and not \{\}

		end %}}}
	end
end
