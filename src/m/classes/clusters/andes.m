%ANDES (Dartmouth) cluster class definition
%
%   Usage:
%      cluster=andes();
%      cluster=andes('np',3);
%      cluster=andes('np',3,'login','username');

classdef andes
	properties (SetAccess=public)
		% {{{
		name          = 'andes'
		login         = '';
		numnodes      = 1;
		cpuspernode   = 16;
		codepath      = '';
		executionpath = '';
		interactive   = 0;
		time          = 10; %in hours
		memory        = 32;  %in Gb
		email         = 'END,FAIL';
		email_domain  = 'dartmouth.edu';
	end
	%}}}
	methods
		function cluster=andes(varargin) % {{{

			%initialize cluster using default settings if provided
			if (exist('andes_settings')==2), andes_settings; end

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
			disp(sprintf('    email: %s (notifications: BEGIN,END,FAIL)',cluster.email));
			disp(sprintf('    codepath:      %s',cluster.codepath));
			disp(sprintf('    executionpath: %s',cluster.executionpath));
			disp(sprintf('    interactive: %i',cluster.interactive));
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
		function BuildKrigingQueueScript(cluster,dirname,modelname,solution,io_gather,isvalgrind,isgprof,isdakota,isoceancoupling) % {{{

			if(isvalgrind), disp('valgrind not supported by cluster, ignoring...'); end
			if(isgprof),    disp('gprof not supported by cluster, ignoring...'); end

			%write queuing script 
			fid=fopen([modelname '.queue'],'w');
			fprintf(fid,'#!/bin/bash\n');
			fprintf(fid,'#SBATCH --job-name=%s\n',modelname);
			fprintf(fid,'#SBATCH --account=ice\n'); %Make sure we use the ICE account for this run
			fprintf(fid,'#SBATCH -o %s.outlog \n',modelname);
			fprintf(fid,'#SBATCH -e %s.errlog \n',modelname);
			fprintf(fid,'#SBATCH --nodes=%i\n',cluster.numnodes);
			fprintf(fid,'#SBATCH --ntasks-per-node=%i\n',cluster.cpuspernode);
			fprintf(fid,'#SBATCH --time=%s\n',datestr(cluster.time/24,'HH:MM:SS')); %walltime is in HH:MM:SS format. cluster.time is in hour
			fprintf(fid,'#SBATCH --mem=%iG\n',cluster.memory);
			if ~isempty(cluster.email)
				fprintf(fid,'#SBATCH --mail-type=%s\n',cluster.email);
				fprintf(fid,'#SBATCH --mail-user=%s@%s\n',cluster.login, cluster.email_domain);
			end
			fprintf(fid,'\n');

			fprintf(fid,'export ISSM_DIR="%s/../"\n',cluster.codepath);
			fprintf(fid,'source $ISSM_DIR/etc/environment.sh\n');      
			fprintf(fid,'cd %s/%s\n\n',cluster.executionpath,dirname);
			fprintf(fid,'srun %s/kriging.exe %s %s\n', cluster.codepath,[cluster.executionpath '/' modelname],modelname);
			if ~io_gather, %concatenate the output files:
				fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
			end
			fclose(fid);
		end
		%}}}
		function BuildQueueScript(cluster,dirname,modelname,solution,io_gather,isvalgrind,isgprof,isdakota,isoceancoupling) % {{{

			if(isvalgrind), disp('valgrind not supported by cluster, ignoring...'); end
			if(isgprof),    disp('gprof not supported by cluster, ignoring...'); end

			%write queuing script
			fid=fopen([modelname '.queue'],'w');
			fprintf(fid,'#!/bin/bash -l\n');
			fprintf(fid,'#SBATCH --job-name=%s\n',modelname);
			fprintf(fid,'#SBATCH --account=ice\n'); %Make sure we use the ICE account for this run
			fprintf(fid,'#SBATCH -o %s.outlog \n',modelname);
			fprintf(fid,'#SBATCH -e %s.errlog \n',modelname);
			fprintf(fid,'#SBATCH --nodes=%i\n',cluster.numnodes);
			fprintf(fid,'#SBATCH --ntasks-per-node=%i\n',cluster.cpuspernode);
			fprintf(fid,'#SBATCH --time=%s\n',eraseBetween(datestr(cluster.time/24,'dd-HH:MM:SS'),1,1)); %walltime is in d-HH:MM:SS format. cluster.time is in hour
			fprintf(fid,'#SBATCH --mem=%iG\n',cluster.memory);
			if ~isempty(cluster.email)
				fprintf(fid,'#SBATCH --mail-type=%s\n',cluster.email);
				fprintf(fid,'#SBATCH --mail-user=%s@%s\n',cluster.login, cluster.email_domain);
			end
			fprintf(fid,'\n');
			fprintf(fid,'mpirun -n %i %s/issm.exe %s %s %s\n',cluster.nprocs(), cluster.codepath,solution,[cluster.executionpath '/' dirname],modelname);
			if ~io_gather, %concatenate the output files:
				fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
			end
			fclose(fid);

			%in interactive mode, create a run file, and errlog and outlog file
			if cluster.interactive,
				fid=fopen([modelname '.run'],'w');
				fprintf(fid,'mpirun -n %i %s/issm.exe %s %s %s\n',cluster.nprocs(), cluster.codepath,solution,[cluster.executionpath '/' dirname],modelname);
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
		function UploadQueueJob(cluster,modelname,dirname,filelist) % {{{

			%compress the files into one zip.
			compressstring=['tar -zcf ' dirname '.tar.gz '];
			for i=1:numel(filelist),
				compressstring = [compressstring ' ' filelist{i}];
			end
			if cluster.interactive,
				compressstring = [compressstring ' ' modelname '.errlog ' modelname '.outlog '];
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
					' && cd ' dirname ' && mv ../' dirname '.tar.gz ./ && tar -zxf ' dirname '.tar.gz  && hostname && sbatch ' modelname '.queue '];
			end
			issmssh(cluster.name,cluster.login,0,launchcommand);
		end %}}}
		function Download(cluster,dirname,filelist) % {{{

			%copy files from cluster to current directory
			directory=[cluster.executionpath '/' dirname '/'];
			issmscpin(cluster.name,cluster.login,0,directory,filelist);

		end %}}}
	end
end
