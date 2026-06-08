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
			port=0; if isprop(cluster,'port'), port=cluster.port; end
			issmscpout(cluster.name,cluster.executionpath,cluster.login,port,{[dirname '.tar.gz']});
		end %}}}
		function LaunchQueueJobSbatch(cluster, modelname, dirname, filelist, restart, batch, format) % {{{
			% on cluster: open tar.gz file and submit job
			%
			% format 1: no scheduler, just source .queue
			% format 2: SLURM, use sbatch
			% format 3: PBS, use qsub

			if format==1
				%No scheduler: source .queue directly
				if ~isempty(restart)
					launchcommand=['source ' cluster.etcpath '/environment.sh && cd ' cluster.executionpath ' && cd ' dirname ' && source ' modelname '.queue '];
				else
					if ~batch
						launchcommand=['source ' cluster.etcpath '/environment.sh && source ' cluster.executionpath '/' dirname '/'  modelname '.queue '];
					else
						launchcommand=['source ' cluster.etcpath '/environment.sh && cd ' cluster.executionpath ' && rm -rf ./' dirname ' && mkdir ' dirname ...
							' && cd ' dirname ' && mv ../' dirname '.tar.gz ./ && tar -zxf ' dirname '.tar.gz '];
					end
				end

			elseif format==2
				%SLURM sbatch launch: unpack the tar and submit the queue script via SSH.
				if ~isempty(restart)
					launchcommand=['cd ' cluster.executionpath ' && cd ' dirname ' && sbatch ' modelname '.queue'];
				else
					launchcommand=['cd ' cluster.executionpath ' && rm -rf ./' dirname ' && mkdir ' dirname ...
						' && cd ' dirname ' && mv ../' dirname '.tar.gz ./ && tar -zxf ' dirname '.tar.gz && sbatch ' modelname '.queue'];
				end
			elseif format==3
				%BPS: use qsub
				if ~isempty(restart)
					launchcommand=['cd ' cluster.executionpath ' && cd ' dirname ' && /PBS/bin/qsub ' modelname '.queue '];
				else
					launchcommand=['cd ' cluster.executionpath ' && rm -rf ./' dirname ' && mkdir ' dirname ...
						' && cd ' dirname ' && mv ../' dirname '.tar.gz ./ && tar -zxf ' dirname '.tar.gz && /PBS/bin/qsub ' modelname '.queue '];
				end

			else
				error('format not supported');
			end
			port=0; if isprop(cluster,'port'), port=cluster.port; end
			issmssh(cluster.name,cluster.login,port,launchcommand);

		end %}}}
		function Download(cluster,dirname,filelist) % {{{

			%Copy output files from the cluster back to the current directory via scp.
			directory=[cluster.executionpath '/' dirname '/'];
			port=0; if isprop(cluster,'port'), port=cluster.port; end
			issmscpin(cluster.name,cluster.login,port,directory,filelist);

		end %}}}
	end
end
