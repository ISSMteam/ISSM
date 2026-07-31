%CLUSTER_DEFAULTS  Default method implementations shared across cluster classes.
%
%   Most remote cluster classes (andes, frontera, pfe, etc.) share
%   identical implementations of UploadQueueJob, LaunchQueueJob, and Download.
%   Rather than duplicating those ~30 lines in every file, each cluster can
%   delegate to the corresponding static method here.
%
%   Clusters whose methods differ from these defaults (e.g. bbftp transfers,
%   interactive path variants, ...) should have their own implementations.

classdef cluster_defaults
	methods (Static)
		function UploadQueueJob(cluster,modelname,dirname,filelist) % {{{

			% Compress the execution directory into a tar.gz and scp it to the cluster.
			% filelist contains full paths; tar with -C so only basenames are stored.
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

			%Copy tar.gz file over
			port=0; if isprop(cluster,'port'); port=cluster.port; end
			issmscpout(cluster.name,cluster.executionpath,cluster.login,port,{[dirname '.tar.gz']});
		end %}}}
		function LaunchQueueJobSbatch(cluster, modelname, dirname, filelist, restart, batch, format) % {{{
			% on cluster: open tar.gz file and submit job
			%
			% format 1: no scheduler, just source .queue
			% format 2: SLURM, use sbatch
			% format 3: PBS, use qsub

			%defaults
			sourceetc_str = ['source ' cluster.etcpath '/environment.sh '];
			untar_str     = [' && cd ' cluster.executionpath ' && rm -rf ./' dirname ' && mkdir ' dirname ' && cd ' dirname ' && mv ../' dirname '.tar.gz ./ && tar -zxf ' dirname '.tar.gz '];

			if ~batch
				if format==1
					%No scheduler: source .queue directly
					submit_str = [' && source ' modelname '.queue '];
				elseif format==2
					%SLURM sbatch launch: unpack the tar and submit the queue script via SSH.
					submit_str = [' && sbatch ' modelname '.queue '];
				elseif format==3
					%BPS: use qsub
					submit_str = [' && /PBS/bin/qsub ' modelname '.queue '];
				else
					error('format not supported');
				end
			else
				submit_str = [' '];
			end

			if strcmp(cluster.name, oshostname())
				untar_str = [' '];
				if ~batch && format==1
					%Special case, do not change directory
					submit_str = [' && source ' cluster.executionpath '/' dirname '/'  modelname '.queue '];
				end
			end

			%Prepare launchcommand
			if ~isempty(restart)
				launchcommand=[sourceetc_str ' && cd ' cluster.executionpath '/' dirname submit_str];
			else
				launchcommand=[sourceetc_str untar_str submit_str];
			end

			%Figure out port
			port=0;
			if isprop(cluster,'port')
				port=cluster.port; 
			end

			%execute launch command
			issmssh(cluster.name, cluster.login, port, launchcommand);

		end %}}}
		function Download(cluster,dirname,filelist) % {{{

			%Copy output files from the cluster back to the current directory via scp.
			directory=[cluster.executionpath '/' dirname '/'];
			port=0; if isprop(cluster,'port'), port=cluster.port; end
			issmscpin(cluster.name,cluster.login,port,directory,filelist);

		end %}}}
	end
end
