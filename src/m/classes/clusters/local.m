%LOCAL cluster class definition
%
%   Usage:
%      cluster=local('name','astrid','np',3);
%      cluster=local('name',oshostname(),'np',3,'login','username');

classdef local
	properties (SetAccess=public)
		% {{{
		name = '';
		np            = 1;
		codepath      = [IssmConfig('ISSM_PREFIX') '/bin'];
		etcpath       = [issmdir() '/etc'];
		executionpath = [issmdir() '/execution'];
		verbose       = 1;
		shell         = '/bin/sh';
		%}}}
	end
	methods
		function cluster=local(varargin) % {{{

			%use provided options to change fields
			options=pairoptions(varargin{:});

			%get name
			cluster.name=getfieldvalue(options,'name',oshostname());

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
			disp(sprintf('    etcpath: %s',cluster.etcpath));
			disp(sprintf('    executionpath: %s',cluster.executionpath));
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
			% Which executable are we calling?
			executable='issm.exe'; % default

			if isdakota,
				executable='issm_dakota.exe';
			end

			fid=fopen([modelname '.queue'],'w');
			fprintf(fid,'#!%s\n',cluster.shell);
			fprintf(fid,'mpiexec -np %i %s/%s %s %s %s \n',cluster.np,cluster.codepath,executable,solution,'./',modelname);
			fclose(fid);
		end
		%}}}
		function UploadQueueJob(cluster,modelname,dirname,filelist)% {{{
			%do nothing really.
		end %}}}
		function LaunchQueueJob(cluster,modelname,dirname,filelist,restart,batch)% {{{
			system(['source ' modelname '.queue']);

		end %}}}
		function Download(cluster,dirname,filelist)% {{{
		end %}}}
	end
end
