%ACI class definition
%
%   Usage:
%      cluster=aci();
%      cluster=aci('np',3);
%      cluster=aci('np',3,'login','username');

classdef aci
	properties (SetAccess=public)  
		% {{{
		name='';
		login='';
		nodes=0;
		ppn=0;
		time=0;
		port=0;
		queue='';
		codepath='';
		executionpath='';
		modules={};
	end
	%}}}
	methods
		function cluster=aci(varargin) % {{{

			%initialize cluster using default settings if provided
			if (exist('aci_settings')==2), aci_settings; end

			%use provided options to change fields
			cluster=AssignObjectFields(pairoptions(varargin{:}),cluster);
		end
		%}}}
		function disp(cluster) % {{{
			%  display the object
			disp(sprintf('class ''%s'' object ''%s'' = ',class(cluster),inputname(1)));
			disp(sprintf('    name: %s',cluster.name));
			disp(sprintf('    login: %s',cluster.login));
			disp(sprintf('    port: %i',cluster.port));
			disp(sprintf('    time: %i (in minutes)',cluster.time));
			disp(sprintf('    nodes: %i (number of nodes)',cluster.nodes));
			disp(sprintf('    ppn: %i (process per nodes)',cluster.ppn));
			disp(sprintf('    queue: %s',cluster.queue));
			disp(sprintf('    codepath: %s',cluster.codepath));
			disp(sprintf('    executionpath: %s',cluster.executionpath));
		end
		%}}}
		function md = checkconsistency(cluster,md,solution,analyses) % {{{

			available_queues={'brp106_a_g_sc_default','open'};
			queue_requirements_time=[Inf 48*60];
			queue_requirements_np=[260 260];

			QueueRequirements(available_queues,queue_requirements_time,queue_requirements_np,cluster.queue,cluster.ppn*cluster.nodes,1)

			%Miscellaneous
			if cluster.ppn>20, md = checkmessage(md,'cannot request more that 20 cores per node'); end
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
			fprintf(fid,'#PBS -A %s\n', cluster.queue); %open or brp106....
			fprintf(fid,'#PBS -l nodes=%i:ppn=%i:stmem\n',cluster.nodes,cluster.ppn);
			fprintf(fid,'#PBS -l walltime=%i\n',cluster.time*60); % walltime is in seconds
			fprintf(fid,'#PBS -o %s.outlog\n',(modelname));
			fprintf(fid,'#PBS -e %s.errlog\n',(modelname));
			fprintf(fid,'source ~/.bashrc\n');
			fprintf(fid,'source $ISSM_DIR/etc/environment.sh\n');
			for i=1:numel(cluster.modules)
				fprintf(fid,'module load %s\n', cluster.modules{i});
			end
			fprintf(fid,'mpirun -np %i %s/issm.exe %s %s %s\n',cluster.nodes*cluster.ppn,cluster.codepath,solution,[cluster.executionpath '/' dirname],modelname);
			if ~io_gather, %concatenate the output files:
				fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
			end
			fclose(fid);
		end %}}}
		function UploadQueueJob(cluster,modelname,dirname,filelist)% {{{

			%compress the files into one zip.
			compressstring=['tar -zcf ' dirname '.tar.gz '];
			for i=1:numel(filelist),
				compressstring = [compressstring ' ' filelist{i}];
			end
			system(compressstring);

			%upload input files
			issmscpout(cluster.name,cluster.executionpath,cluster.login,cluster.port,{[dirname '.tar.gz']});

		end %}}}
		function LaunchQueueJob(cluster,modelname,dirname,filelist,restart,batch)% {{{

			%Execute Queue job
			if ~isempty(restart)
				launchcommand=['cd ' cluster.executionpath ' && cd ' dirname ' && hostname && qsub ' modelname '.queue '];
			else
				launchcommand=['cd ' cluster.executionpath ' && rm -rf ./' dirname ' && mkdir ' dirname ...
					' && cd ' dirname ' && mv ../' dirname '.tar.gz ./ && tar -zxf ' dirname '.tar.gz  && hostname && qsub ' modelname '.queue '];
			end
			issmssh(cluster.name,cluster.login,cluster.port,launchcommand);
		end %}}}
		function Download(cluster,dirname,filelist)% {{{

			%copy files from cluster to current directory
			directory=[cluster.executionpath '/' dirname '/'];
			issmscpin(cluster.name,cluster.login,cluster.port,directory,filelist);

		end %}}}
	end
end
