%CCR cluster class definition
%
%   Usage:
%      cluster=ccr();
%      cluster=ccr('np',3);
%      cluster=ccr('np',3,'login','username');

classdef ccr
        properties (SetAccess=public)
                % {{{
                name           = 'vortex.ccr.buffalo.edu';
                login          = '';
                port           = 0;
                cluster        = 'ub-hpc';% or faculty
                partition      = 'general-compute';
                qos            = 'general-compute';
                account        = '';
                time           = 1; %hr
                numnodes       = 1; %number of nodes
                ntasks         = 1; %number of tasks per node
                cpuspertask    = 1;%number of cpus per task
        exclusive      = false;
                mem            = 1*1000; %GB
                jobname        = '';
                modules        = {'ccrsoft/2023.01' 'cmake/3.22.1' 'matlab/2023b' 'gcc/11.2.0'};
                srcpath        = '';
                codepath       = '';
                executionpath  = '';
                interactive    = 0;
                bbftp          = 0;
        email          = '';
        end
        %}}}
        methods
                function cluster=ccr(varargin) % {{{

                        %initialize cluster using default settings if provided
                        if (exist('ccr_settings')==2), ccr_settings; end

                        %use provided options to change fields
                        cluster=AssignObjectFields(pairoptions(varargin{:}),cluster);
                end
                %}}}
                function disp(cluster) % {{{
                        % display the object
                        disp(sprintf('class ''%s'' object ''%s'' = ',class(cluster),inputname(1)));
                        disp(sprintf('    name: %s',cluster.name));
                        disp(sprintf('    login: %s',cluster.login));
                        disp(sprintf('    port: %i',cluster.port));
                        disp(sprintf('    cluster: %s',cluster.cluster));
                        disp(sprintf('    partition: %s', cluster.partition));
                        disp(sprintf('    queue: %s',cluster.qos));
                        disp(sprintf('    account: %s',cluster.account));
                        disp(sprintf('    time: %i',cluster.time));
                        disp(sprintf('    numnodes: %i',cluster.numnodes));
                        disp(sprintf('    ntasks: %i',cluster.ntasks));
                        disp(sprintf('    cpuspertask: %i',cluster.cpuspertask));
                        disp(sprintf('    memory: %i',cluster.mem));
                        disp(sprintf('    jobname: %s',cluster.jobname));
                        disp(sprintf('    modules: %s',strjoin(cluster.modules,', ')));
                        disp(sprintf('    srcpath: %s',cluster.srcpath));
                        disp(sprintf('    codepath: %s',cluster.codepath));
                        disp(sprintf('    executionpath: %s',cluster.executionpath));
                        disp(sprintf('    interactive: %i',cluster.interactive));
                        disp(sprintf('    bbftp: %s',cluster.bbftp));
                end
                %}}}
                function numprocs=nprocs(cluster) % {{{
                        %compute number of processors
                        numprocs=cluster.numnodes*cluster.cpuspertask;
                end
                %}}}
                function md = checkconsistency(cluster,md,solution,analyses) % {{{
                        if strcmpi(cluster.cluster,'ub-hpc')
                                if ~strcmpi(cluster.qos,cluster.partition)
                                        if ~ismember(cluster.qos,{'supporters','mri','nih'})
						md = md.checkmessage('Value of qos should either match value of partition or be set to "supporters", "mri", or "nih"');
                                        end
                                end

                                available_queues={'debug','general-compute','industry','scavenger','viz','faculty'};
                                queue_requirements_time=[1*3600 72*3600 6*3600 6*3600 24*3600 72*3600];
                                queue_requirements_np=[64 64 56 64 64 64];

                            QueueRequirements(available_queues,queue_requirements_time,queue_requirements_np,cluster.partition,cluster.nprocs(),cluster.time)
            elseif strcmpi(cluster.cluster,'faculty')
                                if strcmpi(cluster.account,'')
                                        md = md.checkmessage('please supply valid account when using the faculty cluster');
                                end
                                if strcmpi(cluster.partition,'sophien')==0 & strcmpi(cluster.qos,'sophien')==0 & strcmpi(cluster.account,sophien')==0
                                        md = md.checkmessage('combination of partition and qos and account invalid');
                                end
                        else
                                md = md.checkmessage('invalid value for cluster');
                        end

                        %Miscellaneous
                        if isempty(cluster.srcpath), md = checkmessage(md,'srcpath empty'); end
                        if isempty(cluster.codepath), md = checkmessage(md,'codepath empty'); end
                        if isempty(cluster.executionpath), md = checkmessage(md,'executionpath empty'); end

                end
                %}}}
                function BuildQueueScript(cluster, md, filename, executable) % {{{

                        %Get variables from md
                        dirname    = md.private.runtimename;
                        modelname  = md.miscellaneous.name;
                        solution   = md.private.solution;
                        io_gather  = md.settings.io_gather;

                        %checks
                        if(md.debug.gprof); disp('gprof not supported by cluster, ignoring...');
            end

                        %write queuing script
                        fid=fopen(filename, 'w');

                        fprintf(fid,'#!/bin/bash -l\n');
                        fprintf(fid,'#SBATCH --time=%i:00:00\n',cluster.time); %walltime is in seconds now converted to hours
                        fprintf(fid,'#SBATCH --ntasks=%i\n', cluster.ntasks);
                        fprintf(fid,'#SBATCH --cpus-per-task=%i\n',cluster.cpuspertask);
                        if strcmpi(cluster.cluster,'faculty')
                                fprintf('#SBATCH --constraint="[SAPPHIRE-RAPIDS-IB|ICE-LAKE-IB|CASCADE-LAKE-IB|EMERALD-RAPIDS-IB]"\n');
                        end
                        fprintf(fid,'#SBATCH --mem=%iG\n',cluster.mem);
                        fprintf(fid,'#SBATCH --job-name=%s\n',cluster.jobname);
                        fprintf(fid,'#SBATCH --output=%s/%s/%s.outlog\n',cluster.executionpath,dirname,modelname);
                        fprintf(fid,'#SBATCH --error=%s/%s/%s.errlog\n',cluster.executionpath,dirname,modelname);
                        fprintf(fid,'#SBATCH --partition=%s\n',cluster.partition);
            fprintf(fid,'#SBATCH --qos=%s\n',cluster.qos);
                        fprintf(fid,'#SBATCH --cluster=%s\n',cluster.cluster);
                        if strcmpi(cluster.account,'')
                                fprintf(fid,'#SBATCH --account=%s\n',cluster.account);
            end
                        for i=1:numel(cluster.modules), fprintf(fid,['module load ' cluster.modules{i} '\n']); end
                        fprintf(fid,'export PATH="$PATH:."\n\n');
                        fprintf(fid,'export ISSM_DIR="%s"\n',cluster.srcpath); %FIXME
                        fprintf(fid,'source $ISSM_DIR/etc/environment.sh\n');
                        fprintf(fid,'cd %s/%s/\n\n',cluster.executionpath,dirname);

                        fprintf(fid,'which mpiexec\n\n');

                        fprintf(fid,'mpiexec -np %i %s/%s %s %s/%s %s\n',cluster.nprocs(),cluster.codepath,executable,solution,cluster.executionpath,dirname,modelname);

                        if ~io_gather %concatenate the output files:
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
            %Upload input files
            disp('uploading input file and queueing script')
            issmscpout(cluster.name,cluster.executionpath,cluster.login,cluster.port,{[dirname '.tar.gz']})
                end
                %}}}
                function LaunchQueueJob(cluster,modelname,dirname,filelist,restart,batch) % {{{
                        if ~isempty(restart)
                                        launchcommand=['cd ' cluster.executionpath ' && cd ' dirname ' && ./' modelname '.queue'];
                                else
                                        launchcommand=['cd ' cluster.executionpath ' && rm -rf ./' dirname ' && mkdir ' dirname ' && cd ' dirname ' && mv ../' dirname '.tar.gz ./ && tar -zxf ' dirname '.tar.gz && sbatch ' modelname '.queue'];
                end
            issmssh(cluster.name,cluster.login,cluster.port,launchcommand);

                end
                %}}}
                function Download(cluster,dirname,filelist) % {{{
            % cluster_defaults.Download(cluster,dirname,filelist);
            directory = [cluster.executionpath '/' dirname '/'];
                        issmscpin(cluster.name,cluster.login,cluster.port,directory,filelist,'auto');
                end %}}}
        end
end
