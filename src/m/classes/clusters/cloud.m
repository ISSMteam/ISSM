%CLOUD cluster class definition
%
%   Usage:
%      cluster=cloud('name','astrid','np',3);
%      cluster=cloud('name',oshostname(),'np',3,'login','username');

classdef cloud
	properties (SetAccess=public)
		% {{{
		name='';
		login='';
		np=1;
		codepath='';
		executionpath='';
		interactive=0;
		%}}}
	end
	methods
		function cluster=cloud(varargin) % {{{

			%initialize cluster using user settings if provided
			if (exist('cloud_settings')==2),
				cloud_settings;
			end

			%OK get other fields
			cluster=AssignObjectFields(pairoptions(varargin{:}),cluster);
		end
		%}}}
		function disp(cluster) % {{{
			%display the object
			disp(sprintf('class ''%s'' object ''%s'' = ',class(cluster),inputname(1)));
			disp(sprintf('    name: %s',cluster.name));
			disp(sprintf('    login: %s',cluster.login));
			disp(sprintf('    np: %i',cluster.np));
			disp(sprintf('    codepath: %s',cluster.codepath));
			disp(sprintf('    executionpath: %s',cluster.executionpath));
			disp(sprintf('    interactive: %i',cluster.interactive));
		end
		%}}}
		function md = checkconsistency(cluster,md,solution,analyses) % {{{
			if cluster.np<1
				md = checkmessage(md,['number of processors should be at least 1']);
			end
			if isnan(cluster.np),
				md = checkmessage(md,'number of processors should not be NaN!');
			end
		end
		%}}}
		function BuildQueueScript(cluster,dirname,modelname,solution,io_gather,isvalgrind,isgprof,isdakota,isoceancoupling) % {{{

			%write queuing script 
			fid=fopen([modelname '.queue'],'w');
			fprintf(fid,'#!/bin/bash\n');
			fprintf(fid,'source %s%s\n',cluster.codepath,'/../etc/environment.sh');
			fprintf(fid,'cd %s\n',[cluster.executionpath '/' dirname]);
			fprintf(fid,'mpiexec -np %i -f /home/mpich2.hosts %s/issm.exe %s %s/%s %s 2> %s.errlog > /dev/stdout | tee %s.outlog',cluster.np,cluster.codepath,solution,cluster.executionpath,dirname,modelname,modelname,modelname);
		end
		%}}}
		function UploadQueueJob(cluster,modelname,dirname,filelist) % {{{

			%compress the files into one zip.
			compressstring=['tar -zcf ' dirname '.tar.gz '];
			for i=1:numel(filelist),
				compressstring = [compressstring ' ' filelist{i}];
			end
			system(compressstring);

			if isempty(cluster.login),
				error('cloud BuildQueueScript: login should be supplied!');
			end
			%upload input files
			issmstscpout(cluster.name,cluster.executionpath,cluster.login,{[dirname '.tar.gz']});

		end %}}}
		function LaunchQueueJob(cluster,modelname,dirname,filelist,restart) % {{{

			if cluster.interactive, 
				disp('sending files to remote cluster. once done, please log into cluster and launch job');
				if ~isempty(restart)
					launchcommand=['cd ' cluster.executionpath ' && cd ' dirname];
				else
					launchcommand=['cd ' cluster.executionpath ' && rm -rf ./' dirname ' && mkdir ' dirname ...
						' && cd ' dirname ' && mv ../' dirname '.tar.gz ./ && tar -zxf ' dirname '.tar.gz'];
				end
			else
				%Execute Queue job
				if ~isempty(restart)
					launchcommand=['cd ' cluster.executionpath ' && cd ' dirname ' && qsub ' modelname '.queue'];
				else
					launchcommand=['cd ' cluster.executionpath ' && rm -rf ./' dirname ' && mkdir ' dirname ...
						' && cd ' dirname ' && mv ../' dirname '.tar.gz ./ && tar -zxf ' dirname '.tar.gz && qsub ' modelname '.queue'];
				end
			end
			issmstssh(cluster.name,cluster.login,launchcommand);
		end %}}}
		function Download(cluster,dirname,filelist) % {{{

			%copy files from cluster to current directory
			directory=[cluster.executionpath '/' dirname '/'];
			issmstscpin(cluster.name,cluster.login,directory,filelist);
		end %}}}
	end
end
