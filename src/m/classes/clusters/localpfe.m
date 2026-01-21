%LOCALPFE cluster class definition
%
%   Usage:
%      cluster=localpfe('name','astrid','np',3);
%      cluster=localpfe('name',oshostname(),'np',3,'login','username');

classdef localpfe
	properties (SetAccess=public) 
		% {{{
		name          = '';
		login         = '';
		np            = 1;
		npocean       = 0;
		port          = 0;
		interactive   = 1;
		codepath      = [IssmConfig('ISSM_PREFIX') '/bin'];
		etcpath       = [issmdir() '/etc'];
		executionpath = [issmdir() '/execution'];
		valgrind      = [issmdir() '/externalpackages/valgrind/install/bin/valgrind'];
		valgrindlib   = [issmdir() '/externalpackages/valgrind/install/lib/libmpidebug.so'];
		valgrindsup   = [issmdir() '/externalpackages/valgrind/issm.supp'];
		verbose       = 1;
		shell         = '/bin/sh';
		%}}}
	end
	methods
		function cluster=localpfe(varargin) % {{{

			%Change the defaults if ispc and not ismingw
			if ispc & ~ismingw,
				cluster.codepath      = [issmdir() '\bin'];
				cluster.etcpath       = [issmdir() '\etc'];
				cluster.executionpath = [issmdir() '\execution'];
			end

			%use provided options to change fields
			options=pairoptions(varargin{:});

			%get name
			cluster.name=getfieldvalue(options,'name',oshostname());

			%initialize cluster using user settings if provided
			if (exist([cluster.name '_settings'])==2), eval([cluster.name '_settings']); end

			%OK get other fields
			cluster=AssignObjectFields(pairoptions(varargin{:}),cluster);
		end
		%}}}
		function disp(cluster) % {{{
			%  display the object
			disp(sprintf('class ''%s'' object ''%s'' = ',class(cluster),inputname(1)));
			disp(sprintf('    name: %s',cluster.name));
			disp(sprintf('    login: %s',cluster.login));
			disp(sprintf('    np: %i',cluster.np));
			disp(sprintf('    npocean: %i',cluster.npocean));
			disp(sprintf('    port: %i',cluster.port));
			disp(sprintf('    interactive: %i',cluster.interactive));
			disp(sprintf('    codepath: %s',cluster.codepath));
			disp(sprintf('    etcpath: %s',cluster.etcpath));
			disp(sprintf('    executionpath: %s',cluster.executionpath));
			disp(sprintf('    valgrind: %s',cluster.valgrind));
			disp(sprintf('    valgrindlib: %s',cluster.valgrindlib));
			disp(sprintf('    valgrindsup: %s',cluster.valgrindsup));
			disp(sprintf('    verbose: %s',cluster.verbose));
			disp(sprintf('    shell: %s',cluster.shell));
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
			%what is the executable being called? 
			executable='issm.exe';
			if isdakota,
				version=IssmConfig('_DAKOTA_VERSION_'); version=str2num(version(1:3));
				if (version>=6),
					executable='issm_dakota.exe';
				end
			end

			fid=fopen([modelname '.queue'],'w');
			fprintf(fid,'#!%s\n',cluster.shell);
			fprintf(fid,'mpiexec -np %i %s/%s %s %s %s \n',cluster.np,cluster.codepath,executable,solution,cluster.executionpath,modelname);
			fclose(fid);

			%in interactive mode, create a run file, and errlog and outlog file
			if cluster.interactive,
				fid=fopen([modelname '.run'],'w');
				if cluster.interactive==10,
						fprintf(fid,'module unload mpi-mvapich2/1.4.1/gcc\n');
						fprintf(fid,'mpiexec -np %i %s/%s %s %s %s\n',cluster.np,cluster.codepath,executable,solution,[pwd() '/run'],modelname);
				else
					if ~isvalgrind,
						fprintf(fid,'mpiexec -np %i %s/%s %s %s %s\n',cluster.np,cluster.codepath,executable,solution,cluster.executionpath,modelname);
						%fprintf(fid,'mpiexec -np %i %s/%s %s %s %s\n',cluster.nprocs(),cluster.codepath,executable,solution,[cluster.executionpath '/Interactive' num2str(cluster.interactive)],modelname);
					else
						fprintf(fid,'mpiexec -np %i valgrind --leak-check=full %s/%s %s %s %s\n',cluster.np,cluster.codepath,executable,solution,[cluster.executionpath '/Interactive' num2str(cluster.interactive)],modelname);
					end
				end
				if ~io_gather, %concatenate the output files:
					fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
				end
				fclose(fid);

				fid=fopen([modelname '.errlog'],'w'); fclose(fid);
				fid=fopen([modelname '.outlog'],'w'); fclose(fid);
			end
		end
		%}}}
		function BuildQueueScriptMultipleModels(cluster,dirname,modelname,solution,dirnames,modelnames,nps) % {{{
		
			%some checks: 
			if isempty(modelname), error('BuildQueueScriptMultipleModels error message: need a non empty model name!');end

			%what is the executable being called? 
			executable='issm_slc.exe';

			if ispc & ~ismingw, error('BuildQueueScriptMultipleModels not support yet on windows machines');end;
			
			%write queuing script 
			fid=fopen([modelname '.queue'],'w');
			
			fprintf(fid,'#!%s\n',cluster.shell);

			%number of cpus: 
			mpistring=sprintf('mpiexec -np %i ',cluster.np);

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

			%log files: 
			if ~cluster.interactive,
				mpistring=[mpistring sprintf('2> %s.errlog> %s.outlog',modelname,modelname)];
			end

			%write this long string to disk: 
			fprintf(fid,mpistring);
			fclose(fid);

			%in interactive mode, create a run file, and errlog and outlog file
			if cluster.interactive,
				fid=fopen([modelname '.errlog'],'w'); fclose(fid);
				fid=fopen([modelname '.outlog'],'w'); fclose(fid);
			end
		end
		%}}}
		function BuildQueueScriptIceOcean(cluster,dirname,modelname,solution,io_gather,isvalgrind,isgprof,isdakota) % {{{

			%write queuing script 
			%what is the executable being called? 
			executable='issm_ocean.exe';

			fid=fopen([modelname '.queue'],'w');
			fprintf(fid,'#!%s\n',cluster.shell);
			fprintf(fid,'mpiexec -np %i %s/%s %s %s %s : -np %i ./mitgcmuv\n',cluster.np,cluster.codepath,executable,solution,cluster.executionpath,modelname,cluster.npocean);
			fclose(fid);

			%in interactive mode, create a run file, and errlog and outlog file
			if cluster.interactive,
				fid=fopen([modelname '.errlog'],'w'); fclose(fid);
				fid=fopen([modelname '.outlog'],'w'); fclose(fid);
			end
		end
		%}}}
		function BuildKrigingQueueScript(cluster,modelname,solution,io_gather,isvalgrind,isgprof) % {{{

			%write queuing script 
			if ~ispc,

				fid=fopen([modelname '.queue'],'w');
				fprintf(fid,'#!/bin/sh\n');
				if ~isvalgrind,
					if cluster.interactive
						fprintf(fid,'mpiexec -np %i %s/kriging.exe %s %s ',cluster.np,cluster.codepath,[cluster.executionpath '/' modelname],modelname);
					else
						fprintf(fid,'mpiexec -np %i %s/kriging.exe %s %s 2> %s.errlog >%s.outlog ',cluster.np,cluster.codepath,[cluster.executionpath '/' modelname],modelname,modelname,modelname);
					end
				elseif isgprof,
					fprintf(fid,'\n gprof %s/kriging.exe gmon.out > %s.performance',cluster.codepath,modelname);
				else
					%Add --gen-suppressions=all to get suppression lines
					fprintf(fid,'LD_PRELOAD=%s \\\n',cluster.valgrindlib);
					fprintf(fid,'mpiexec -np %i %s --leak-check=full --suppressions=%s %s/kriging.exe %s %s 2> %s.errlog >%s.outlog ',...
						cluster.np,cluster.valgrind,cluster.valgrindsup,cluster.codepath,[cluster.executionpath '/' modelname],modelname,modelname,modelname);
				end
				if ~io_gather, %concatenate the output files:
					fprintf(fid,'\ncat %s.outbin.* > %s.outbin',modelname,modelname);
				end
				fclose(fid);

			else % Windows

				fid=fopen([modelname '.bat'],'w');
				fprintf(fid,'@echo off\n');
				if cluster.interactive
					fprintf(fid,'"%s/issm.exe" %s "%s" %s ',cluster.codepath,solution,[cluster.executionpath '/' modelname],modelname);
				else
					fprintf(fid,'"%s/issm.exe" %s "%s" %s 2> %s.errlog >%s.outlog',...
						cluster.codepath,solution,[cluster.executionpath '/' modelname],modelname,modelname,modelname);
				end
				fclose(fid);
			end

			%in interactive mode, create a run file, and errlog and outlog file
			if cluster.interactive,
				fid=fopen([modelname '.errlog'],'w'); fclose(fid);
				fid=fopen([modelname '.outlog'],'w'); fclose(fid);
			end
		end
		%}}}
		function UploadQueueJob(cluster,modelname,dirname,filelist)% {{{
			if ~ispc | ismingw,

				%compress the files into one zip.
				compressstring=['tar -zcf ' dirname '.tar.gz '];
				for i=1:numel(filelist),
					compressstring = [compressstring ' ' filelist{i}];
				end
				if cluster.interactive,
					compressstring = [compressstring ' ' modelname '.run '  modelname '.errlog ' modelname '.outlog '];
				end
				system(compressstring);

				issmscpout(cluster.name,cluster.executionpath,cluster.login,cluster.port,{[dirname '.tar.gz']});
			end
		end %}}}
		function LaunchQueueJob(cluster,modelname,dirname,filelist,restart,batch)% {{{

			%figure out what shell extension we will use:
			if isempty(strfind(cluster.shell,'csh')),
				shellext='sh';
			else
				shellext='csh';
			end

			if cluster.verbose, %Execute Queue job end

			launchcommand=['cd ' cluster.executionpath ' && rm -rf *.lock && rm -rf ADOLC* && tar -zxf ' dirname '.tar.gz  && rm -rf *.tar.gz'];
			issmssh(cluster.name,cluster.login,cluster.port,launchcommand);

		end %}}}
		function LaunchQueueJobIceOcean(cluster,modelname,dirname,filelist,restart,batch)% {{{

			%figure out what shell extension we will use:
			if isempty(strfind(cluster.shell,'csh')),
				shellext='sh';
			else
				shellext='csh';
			end

			if cluster.verbose, %Execute Queue job end

			launchcommand=['cd ' cluster.executionpath ' && rm -rf *.lock && tar -zxf ' dirname '.tar.gz  && rm -rf *.tar.gz'];
			issmssh(cluster.name,cluster.login,cluster.port,launchcommand);

		end %}}}
		function Download(cluster,dirname,filelist)% {{{

			if ispc && ~ismingw,
				%do nothing
				return;
			end

			%copy files from cluster to current directory
			issmscpin(cluster.name,cluster.login,cluster.port,cluster.executionpath,filelist);
		end %}}}
	end
end
