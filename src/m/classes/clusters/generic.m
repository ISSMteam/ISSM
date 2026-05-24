%GENERIC cluster class definition
%
%   Usage:
%      cluster=generic('name','astrid','np',3);
%      cluster=generic('name',oshostname(),'np',3,'login','username');
%
%   TODO:
%   - Add support for restart to Windows (under MSYS2), then activate tests 125 
%   and 126 in test suite
%

classdef generic
	properties (SetAccess=public)
		% {{{
		name          = '';
		login         = '';
		np            = 1;
		npocean       = 1;
		port          = 0;
		interactive   = 1;
		codepath      = [issmdir() '/bin'];
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
		function cluster=generic(varargin) % {{{

			%Change the defaults if ispc
			if ispc
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
			disp(sprintf('    executionpath: %s',cluster.executionpath));
			disp(sprintf('    etcpath: %s',cluster.etcpath));
			disp(sprintf('    valgrind: %s',cluster.valgrind));
			disp(sprintf('    valgrindlib: %s',cluster.valgrindlib));
			disp(sprintf('    valgrindsup: %s',cluster.valgrindsup));
			disp(sprintf('    verbose: %i',cluster.verbose));
			disp(sprintf('    shell: %s',cluster.shell));
		end
		%}}}
		function numprocs=nprocs(cluster) % {{{
			numprocs=cluster.np;
		end
		%}}}
		function md = checkconsistency(cluster,md,solution,analyses) % {{{
			if cluster.np<1
				md = checkmessage(md,['number of processors should be at least 1']);
			end
			if isnan(cluster.np)
				md = checkmessage(md,'number of processors should not be NaN!');
			end
		end
		%}}}
		function BuildQueueScript(cluster, md, filename, executable) % {{{

			% Unpack fields used below
			dirname         = md.private.runtimename;
			modelname       = md.miscellaneous.name;
			solution        = md.private.solution;
			io_gather       = md.settings.io_gather;
			isvalgrind      = md.debug.valgrind;
			isgprof         = md.debug.gprof;
			isdakota        = md.qmu.isdakota;
			isoceancoupling = md.transient.isoceancoupling;

			% Determine which executable to call
			executable = 'issm.exe';
			if isdakota
				dakota_version_str = IssmConfig('_DAKOTA_VERSION_'); version = str2num(dakota_version_str(1:3));
				if version >= 6, executable = 'issm_dakota.exe'; end
			end
			if isoceancoupling
				executable = 'issm_ocean.exe';
			end

			if ~ispc()

				% Verify the executable exists
				exepath = [cluster.codepath '/' executable];
				if ~exist(exepath, 'file')
					error('File %s does not exist', exepath);
				end

				% Escape spaces in codepath for the shell script
				codepath = strrep(cluster.codepath, ' ', '\ ');
				execpath = [cluster.executionpath '/' dirname];

				% Build the mpi prefix once (empty string when MPI is not available)
				mpiprefix = '';
				if IssmConfig('_HAVE_MPI_')
					mpiprefix = sprintf('mpiexec -np %i ', cluster.np);
				end

				% Build the core command string
				if isvalgrind
					if ismac
						vgflags = '--leak-check=full --show-leak-kinds=all --error-limit=no --dsymutil=yes';
					else
						vgflags = '--leak-check=full --error-limit=no';
					end
					cmd = sprintf('%s%s %s --suppressions=%s %s/%s %s %s %s 2> %s.errlog > %s.outlog', ...
						mpiprefix, cluster.valgrind, vgflags, cluster.valgrindsup, ...
						codepath, executable, solution, execpath, modelname, modelname, modelname);

				elseif isgprof
					cmd = sprintf('gprof %s/issm.exe gmon.out > %s.performance', cluster.codepath, modelname);

				elseif cluster.interactive
					cmd = sprintf('%s%s/%s %s %s %s', mpiprefix, codepath, executable, solution, execpath, modelname);

				else
					% Non-interactive: redirect output and run in background
					if IssmConfig('_HAVE_MPI_')
						cmd = sprintf('%s%s/%s %s %s %s 2> %s/%s.errlog > %s/%s.outlog &', ...
							mpiprefix, codepath, executable, solution, execpath, modelname, execpath, modelname, execpath, modelname);
					else
						cmd = sprintf('%s/%s %s %s %s 2> %s.errlog > %s.outlog &', ...
							codepath, executable, solution, execpath, modelname, modelname, modelname);
					end
				end

				% Write the queue script
				fid = fopen(filename, 'w');
				fprintf(fid, '#!%s\n', cluster.shell);
				fprintf(fid, '%s', cmd);
				if ~io_gather
					% Concatenate distributed output files into one
					fprintf(fid, '\ncat %s.outbin.* > %s.outbin', modelname, modelname);
				end
				fclose(fid);

			else % Windows

				batfilename = [filename(1:end-6) '.bat'];
				execdir     = [cluster.executionpath '\' dirname];
				fid         = fopen(batfilename, 'w');
				fprintf(fid, '@echo off\n');
				if cluster.np > 1
					fprintf(fid, '"C:\\Program Files\\Microsoft MPI\\Bin\\mpiexec.exe" -n %i "%s\\%s" %s "%s" %s', ...
						cluster.np, cluster.codepath, executable, solution, execdir, modelname);
				else
					fprintf(fid, '"%s\\%s" %s "%s" %s', cluster.codepath, executable, solution, execdir, modelname);
				end
				fclose(fid);

			end
		end
		%}}}
		function BuildQueueScriptMultipleModels(cluster, slm, dirnames, modelnames, nps, filename) % {{{

         %Get variables from slm
         dirname         = slm.private.runtimename;
         modelname       = slm.miscellaneous.name;
         solution        = slm.private.solution;

         %checks
			if(isempty(modelname)) error('BuildQueueScriptMultipleModels error message: need a non empty model name!');end
			if(ispc()) error('BuildQueueScriptMultipleModels not support yet on windows machines');end;

			%what is the executable being called?
			executable='issm_slc.exe';

			%write queuing script
			fid=fopen(filename, 'w');

			fprintf(fid,'#!%s\n',cluster.shell);

			%number of cpus:
			mpistring = sprintf('mpiexec -np %i ',cluster.np);

			%executable:
			mpistring=[mpistring sprintf('%s/%s ',cluster.codepath,executable)];

			%solution name:
			mpistring=[mpistring sprintf('%s ',solution)];

			%execution directory and model name:
			mpistring=[mpistring sprintf('%s/%s %s',cluster.executionpath,dirname,modelname)];

			%inform main executable of how many icecaps, glaciers and earth models are being run:
			mpistring=[mpistring sprintf(' %i ',length(dirnames))];

			%icecaps, glaciers and earth location, names and number of processors associated:
			for i=1:length(dirnames)
			mpistring=[mpistring sprintf(' %s/%s %s %i ',cluster.executionpath,dirnames{i},modelnames{i},nps{i})];
			end

			%log files:
			if ~cluster.interactive
				mpistring=[mpistring sprintf('2> %s.errlog> %s.outlog',modelname,modelname)];
			end

			%write this long string to disk:
			fprintf(fid,mpistring);
			fclose(fid);
		end
		%}}}
		function BuildQueueScriptIceOcean(cluster, md, filename) % {{{

         %Get variables from md
         dirname         = md.private.runtimename;
         modelname       = md.miscellaneous.name;
         solution        = md.private.solution;
         io_gather       = md.settings.io_gather;
         isvalgrind      = md.debug.valgrind;
         isgprof         = md.debug.gprof;
         isdakota        = md.qmu.isdakota;
         isoceancoupling = md.transient.isoceancoupling;

			%what is the executable being called?
			executable='issm_ocean.exe';

			fid=fopen(filename, 'w');
			fprintf(fid,'#!%s\n',cluster.shell);
			if ~isvalgrind
				fprintf(fid,'mpiexec -np %i %s/%s %s %s %s : -np %i ./mitgcmuv\n',cluster.np,cluster.codepath,executable,solution,cluster.executionpath,modelname,cluster.npocean);

			else
				fprintf(fid,'mpiexec -np %i %s --leak-check=full --error-limit=no --dsymutil=yes --suppressions=%s  %s/%s %s %s %s : -np %i ./mitgcmuv\n',...
					cluster.np,cluster.valgrind,cluster.valgrindsup,cluster.codepath,executable,solution,cluster.executionpath,modelname,cluster.npocean);
			end
			fclose(fid);
		end
		%}}}
		function BuildKrigingQueueScript(cluster, md, filename) % {{{

         %Get variables from md
         dirname         = md.private.runtimename;
         modelname       = md.miscellaneous.name;
         solution        = md.private.solution;
         io_gather       = md.settings.io_gather;
         isvalgrind      = md.debug.valgrind;
         isgprof         = md.debug.gprof;
         isdakota        = md.qmu.isdakota;
         isoceancoupling = md.transient.isoceancoupling;

			%write queuing script
			if ~ispc()

				fid=fopen(filename, 'w');
				fprintf(fid,'#!/bin/sh\n');
				if ~isvalgrind
					if cluster.interactive
						fprintf(fid,'mpiexec -np %i %s/kriging.exe %s %s ',cluster.np,cluster.codepath,[cluster.executionpath '/' modelname],modelname);
					else
						fprintf(fid,'mpiexec -np %i %s/kriging.exe %s %s 2> %s.errlog >%s.outlog ',cluster.np,cluster.codepath,[cluster.executionpath '/' modelname],modelname,modelname,modelname);
					end
				elseif isgprof
					fprintf(fid,'\n gprof %s/kriging.exe gmon.out > %s.performance',cluster.codepath,modelname);
				else
					%Add --gen-suppressions=all to get suppression lines
					%fprintf(fid,'LD_PRELOAD=%s \\\n',cluster.valgrindlib); it could be deleted
					fprintf(fid,'mpiexec -np %i %s --leak-check=full --suppressions=%s %s/kriging.exe %s %s 2> %s.errlog >%s.outlog ',...
						cluster.np,cluster.valgrind,cluster.valgrindsup,cluster.codepath,[cluster.executionpath '/' modelname],modelname,modelname,modelname);
				end
				fclose(fid);
			else % Windows
				error('not supported');
			end
		end
		%}}}
		function UploadQueueJob(cluster,modelname,dirname,filelist) % {{{

			if ~ispc

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

				issmscpout(cluster.name,cluster.executionpath,cluster.login,cluster.port,{[dirname '.tar.gz']});
			end
		end %}}}
		function LaunchQueueJob(cluster,modelname,dirname,filelist,restart,batch) % {{{

			if ~ispc
				cluster_defaults.LaunchQueueJobSbatch(cluster,modelname,dirname,filelist,restart,batch, 1);
			else
				batfile=[cluster.executionpath '\' dirname '\' modelname '.bat'];
				system(['"' batfile '"']);
			end

		end %}}}
		function Download(cluster,dirname,filelist) % {{{

			%copy files from cluster to current directory
			directory=[cluster.executionpath '/' dirname '/'];
			issmscpin(cluster.name,cluster.login,cluster.port,directory,filelist);
		end %}}}
	end
end
