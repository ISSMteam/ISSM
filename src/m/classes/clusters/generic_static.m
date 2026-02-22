%GENERIC cluster class definition
%
%   Usage:
%      cluster=generic_static('name','astrid','np',3);
%      cluster=generic('name',oshostname(),'np',3,'login','username');
%
%   TODO:
%   - Add support for restart to Windows (under MSYS2), then activate tests 125 
%   and 126 in test suite
%

classdef generic_static
	properties (SetAccess=public)
		% {{{
		name='';
		np=1;
		codepath=fileparts(which('issm.exe'));
		executionpath='.';
		interactive=1;
		shell='/bin/sh';
		%}}}
	end
	methods
		function cluster=generic_static(varargin) % {{{

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
			disp(sprintf('    np: %i',cluster.np));
			disp(sprintf('    codepath: %s',cluster.codepath));
			disp(sprintf('    executionpath: %s',cluster.executionpath));
			disp(sprintf('    interactive: %s',cluster.interactive));
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
			if isnan(cluster.np),
				md = checkmessage(md,'number of processors should not be NaN!');
			end
		end
		%}}}
		function BuildQueueScript(cluster,dirname,modelname,solution,io_gather,isvalgrind,isgprof,isdakota,isoceancoupling) % {{{
			% Which executable are we calling?
			executable='issm.exe'; % default

			if isdakota,
				version=IssmConfig('_DAKOTA_VERSION_');
				version=str2num(version(1:3));
				if (version>=6),
					executable='issm_dakota.exe';
				end
			end
			if isoceancoupling,
				executable='issm_ocean.exe';
			end

			if ~ispc(),
				% Check that executable exists at the right path
				if ~exist([cluster.codepath '/' executable],'file'),
					error(['File ' cluster.codepath '/' executable ' does not exist']);
				end

				% Process codepath and prepend empty spaces with \ to avoid errors in queuing script
				codepath=strrep(cluster.codepath,' ','\ ');

				% Write queuing script
				fid=fopen([modelname '.queue'],'w');
				fprintf(fid,'#!%s\n',cluster.shell);
				fprintf(fid,['%s/mpiexec -np %i %s/%s %s %s %s \n'],codepath,cluster.np,codepath,executable,solution,'./',modelname);
				fclose(fid);
			else % Windows
				fid=fopen([modelname '.bat'],'w');
				fprintf(fid,'@echo off\n');

				if cluster.np>1,
					fprintf(fid,'"%s\\mpiexec.exe" -n %i "%s/%s" %s ./ %s',cluster.codepath,cluster.np,cluster.codepath,executable,solution,modelname);
				else
					fprintf(fid,'"%s\\%s" %s ./ %s',cluster.codepath,executable,solution,modelname);
				end
				fclose(fid);
			end

			%Create an errlog and outlog file
			fid=fopen([modelname '.errlog'],'w');
			fclose(fid);
			fid=fopen([modelname '.outlog'],'w');
			fclose(fid);
		end
		%}}}
		function UploadQueueJob(cluster,modelname,dirname,filelist) % {{{
			% Do nothing
			return;
		end %}}}
		function LaunchQueueJob(cluster,modelname,dirname,filelist,restart,batch) % {{{
			if ~ispc,
				% Figure out which file extension to use
				if isempty(strfind(cluster.shell,'csh')),
					shellext='sh';
				else
					shellext='csh';
				end

				launchcommand=['source  ' modelname '.queue '];
				issmssh(cluster.name,'',0,launchcommand);
			else
				system([modelname '.bat']);
			end
		end %}}}
		function Download(cluster,dirname,filelist) % {{{
			% Do nothing
			return;
		end %}}}
	end
end
