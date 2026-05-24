%COMPUTECANADA class definition
%
%   Usage:
%      cluster=computecanada();
%      cluster=computecanada('np',3);
%      cluster=computecanada('np',3,'login','username');

classdef computecanada
    properties (SetAccess=public)  
		 % {{{
		 name           = ''
		 login          = '';
		 numtasks       = 1;
		 cpuspertask    = 8;
		 port           = 0;
		 projectaccount = '';
		 codepath       = '';
		 executionpath  = '';
		 time           = 24*60;
		 memory         = 2;
		 email          = '';
		 mailtype       = '';
	 end
	 %}}}
	 methods
		 function cluster=computecanada(varargin) % {{{

			 %initialize cluster using default settings if provided
			 if (exist('computecanada_settings')==2), computecanada_settings; end

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
			 disp(sprintf('    numtasks: %i',cluster.numtasks));
			 disp(sprintf('    cpuspertask: %i',cluster.cpuspertask));
			 disp(sprintf('    projectaccount: %s',cluster.projectaccount));
			 disp(sprintf('    codepath: %s',cluster.codepath));
			 disp(sprintf('    executionpath: %s',cluster.executionpath));
			 disp(sprintf('    time: %i',cluster.time));
			 disp(sprintf('    memory: %i',cluster.memory));
			 disp(sprintf('    email: %s', cluster.email));
			 disp(sprintf('    mailtype: %s', cluster.mailtype));
			 
		 end
		 %}}}
		 function numprocs=np(cluster) % {{{
			 %compute number of processors
			 numprocs=cluster.numtasks*cluster.cpuspertask;
		 end
		 %}}}
		 function md = checkconsistency(cluster,md,solution,analyses) % {{{
			 if isempty(cluster.name), md = checkmessage(md,'name empty'); end
			 if isempty(cluster.login), md = checkmessage(md,'login empty'); end
			 if ~(cluster.numtasks > 0), md = checkmessage(md,'numtasks must be > 0'); end
			 if ~(cluster.cpuspertask > 0), md = checkmessage(md,'cpuspertask must be > 0'); end
			 if ~(cluster.port >= 0), md = checkmessage(md,'port must be >= 0'); end
			 if isempty(cluster.email), md = checkmessage(md,'email empty'); end
			 if isempty(cluster.mailtype), md = checkmessage(md,'mailtype empty'); end
			 if isempty(cluster.projectaccount), md = checkmessage(md,'projectaccount empty'); end
			 if isempty(cluster.codepath), md = checkmessage(md,'codepath empty'); end
			 if isempty(cluster.executionpath), md = checkmessage(md,'executionpath empty'); end
			 if ~(cluster.time > 0), md = checkmessage(md,'time must be > 0'); end
			 if ~(cluster.memory > 0), md = checkmessage(md,'memory must be > 0'); end
		 end
		 %}}}
		 function BuildQueueScript(cluster, md, filename, executable) % {{{

			%Get variables from md
			dirname   = md.private.runtimename;
			modelname = md.miscellaneous.name;
			solution  = md.private.solution;
			io_gather = md.settings.io_gather;

			%checks
			 if(md.debug.valgrind) disp('valgrind not supported by cluster, ignoring...'); end
			 if(md.debug.gprof)    disp('gprof not supported by cluster, ignoring...'); end

			 %write queuing script 
			 fid=fopen(filename, 'w');
			 fprintf(fid,'#!/bin/bash\n');
			 fprintf(fid,'#SBATCH --job-name=%s\n',modelname);
			 fprintf(fid,'#SBATCH --account=%s \n',cluster.projectaccount);
			 fprintf(fid,'#SBATCH --ntasks=%i  \n',cluster.numtasks);
			 fprintf(fid,'#SBATCH --cpus-per-task=%i\n',cluster.cpuspertask);
			 fprintf(fid,'#SBATCH --time=%i\n',cluster.time); %walltime is in minutes
			 fprintf(fid,'#SBATCH --mem-per-cpu=%igb\n',cluster.memory); %memory in in gigabytes
			 fprintf(fid,'#SBATCH --mail-user=%s\n',cluster.email); %email
			 fprintf(fid,'#SBATCH --mail-type=%s\n',cluster.mailtype); 
			 fprintf(fid,'#SBATCH --output=%s.outlog \n',modelname);
			 fprintf(fid,'#SBATCH --error=%s.errlog \n\n',modelname);
			 fprintf(fid,'export ISSM_DIR="%s/../"\n',cluster.codepath); %FIXME
			 fprintf(fid,'cd %s/%s\n\n',cluster.executionpath,dirname);
			 fprintf(fid,'srun -n %i %s/%s %s %s %s\n',cluster.np(),cluster.codepath,executable,solution,[cluster.executionpath '/' dirname],modelname); 
			 if ~io_gather, %concatenate the output files:
				 fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
			 fclose(fid);
		 end %}}}
		 function UploadQueueJob(cluster,modelname,dirname,filelist)% {{{
			 cluster_defaults.UploadQueueJob(cluster,modelname,dirname,filelist);
		 end %}}}
		 function LaunchQueueJob(cluster,modelname,dirname,filelist,restart,batch)% {{{
			 cluster_defaults.LaunchQueueJobSbatch(cluster,modelname,dirname,filelist,restart,batch, 2);
		 end %}}}
		 function Download(cluster,dirname,filelist)% {{{
			 cluster_defaults.Download(cluster,dirname,filelist);
		 end %}}}
	end
end
