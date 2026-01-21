%MODELLIST class definition
%
%   Usage:
%      modellist=modellist({md1 md2 md3});

classdef modellist
	properties (SetAccess=public) 
		models  = cell(0,1);
		cluster = generic();
	end
	methods
		function md_list=modelsextract(md,flags,minel,varargin) % {{{
			%modelsextract - extract several self contained models according to a list of element flags.
			%
			%   The difference between this routine and the modelextract.m routine (without an 's') is that 
			%   as many models are extracted as there are closed contours defined in area. 
			%   This routine is needed for example when doing data assimilation of ice shelves in Antarctica. 
			%   Many independent ice shelves are present, and we don't want data assimilation on one ice shelf 
			%   to be hindered by another totally independent ice shelf.
			%
			%   Usage:
			%      md_list=modelsextract(md,elementfalgs,minel);
			%
			%   Examples:
			%      md_list=modelsextract(md...,,1000);
			%
			%   See also: EXTRUDE, COLLAPSE, MODELEXTRACT

			disp('selecting pools of elements');
			%go through flags and build as many independent element flags as there are groups of connected 1s
			%in flags.

			%2D or 3D?
			if dimension(md.mesh)==3,
				numberofelements=md.mesh.numberofelements2d; %this will be forgotten when we get out.
				flags=project2d(md,flags,1);
			else
				numberofelements=md.mesh.numberofelements;
			end

			%recover extra arguments: 
			distance=0;
			if nargin==4,
				distance=varargin{1};
			end

			flag_list=cell(0,1);

			for i=1:size(flags,1),

				if (flags(i)),

					%ok, we are sure element i is part of a new pool.
					pool=zeros(numberofelements,1);
					pool=PropagateFlagsFromConnectivity(md.mesh.elementconnectivity,pool,i,flags);
					flag_list{end+1,1}=pool;

					%speed up rest of computation by taking pool out of flags: 
					pos=find(pool);flags(pos)=0;

				end
			end

			%go through flag_list and discard any pool of less than minel elements: 
			ex_pos=[];
			for i=1:length(flag_list),
				if length(find(flag_list{i}))<minel,
					ex_pos=[ex_pos; i];
				end
			end
			flag_list(ex_pos)=[];

			%now, if distance was specified, expand the flag_list by distance km: 
			if distance,
				for i=1:length(flag_list),
					flag_list{i}=PropagateFlagsUntilDistance(md,flag_list{i},distance);
				end
			end

			%now, go use the pools of flags to extract models: 
			disp(['extracting ' num2str(size(flag_list,1)) ' models']);
			models=cell(0,1);

			for i=1:size(flag_list,1),
				disp(['   ' num2str(i) '/' num2str(size(flag_list,1))]);
				if dimension(md.mesh)==3,
					flags2d=flag_list{i};
					realflags=project3d(md,flags2d,'element');
				else
					realflags=flag_list{i};
				end
				models{end+1,1}=modelextract(md,realflags);
			end

			%return model list
			md_list=modellist(models);

		end %end of this function }}}
		function md_list=modelsextractfromdomains(md,directory) % {{{
			%modelsextractfromdomains- extract several self contained models according to a list of domains
			%
			%   Usage:
			%      md_list=modelsextractfromdomains(md,'Basins/');
			%
			%   Examples:
			%      md_list=modelsextract(md,'Basins/');
			%
			%   See also: MODELSEXTRACTS, MODELEXTRACT

			%go into directory and get list of files.
			cd(directory);
			basins=listfiles;
			cd ..

			models=cell(0,1);
			for i=1:length(basins),
				models{end+1,1}=modelextract(md,[directory '/' basins{i}]);
			end

			%return model list: 
			md_list=modellist(models);

		end % }}}
		function self = modellist(varargin) % {{{

			%initialize list
			if nargin==0,
				%Do nothing,
			elseif nargin==1,
				if ~isa(varargin{1},'cell'),
					error('not supported yet');
				end

				celllist=varargin{1};

				%check on size of cell list: 
				if (size(celllist,2)~=1),
					error('modellist constructor error message: list of models should be a cell list of column size 1');
				end

				%check that only models are in the celllist: 
				for i=1:size(celllist,1),
					if ~isa(celllist{i},'model')
						error(['modellist constructor error message: element ' num2str(i) ' of cell list is not a model!']);
					end
				end

				self.models  = celllist;
				self.cluster = self.models{1}.cluster;
			end
		end % }}}
		function val = get(self, propName)% {{{
		%GET - gets model propertie from a specified object ans returns the value
		% 
		%   Usage:
		%      val = get(a, propName)

			switch propName
				case 'numberofelements'
					val = self.numberofelements;
				case 'numberofnodes'
					val = self.numberofnodes;
				case 'elements' 
					val = self.elements;
				case 'x' 
					val = self.x;
				case 'y' 
					val = self.y;
				case 'z' 
					val = self.z;
				otherwise
					error(['get error message: ' propName,' is not a valid model property'])
			end
		end % }}}
		function self = loadmultipleresultsfromcluster(self) % {{{
			%LOADMULTIPLERESULTSFROMCLUSTER - load multiple results of solution sequences from cluster
			%
			%   Usage:
			%      self=loadresultsfromcluster(self);

			nummodels=length(self.models);

			%Get cluster settings
			cluster=self.cluster;
			name=self.name;
			cluster_rc_location=which('cluster.rc');
			[codepath,executionpath]=ClusterParameters(cluster,cluster_rc_location);

			%Remote tar: 
			disp('tarring results');
			issmssh(cluster,['"cd ' executionpath '/' name ' && rm -rf file_list.txt ModelResults.tar.gz && find -iname ''*-*vs*.outbin'' > file_list.txt && tar zcvf ModelResults.tar.gz --files-from file_list.txt  && rm -rf file_list.txt "']);

			%copy results from cluster to present directory
			scpin(cluster, [executionpath '/' name], {'ModelResults.tar.gz'});

			%untar:
			!tar -zxvf ModelResults.tar.gz

			%ok, go through list and load results from disk: 
			for i=1:nummodels,
				%load  results for this model
				self.models{i}=loadresultsfromdisk(self.models{i},[name '-' num2str(i) 'vs' num2str(nummodels) '.outbin']);

				delete([name '-' num2str(i) 'vs' num2str(nummodels) '.outbin']);
			end

			%erase files 
			delete('ModelResults.tar.gz');
		end % }}}
		function self = solve(self,varargin)% {{{
			%SOLVE - apply solution sequence for  a list of models. Used in batch mode.
			%
			%   Usage:
			%      self=solve(self,varargin)
			%      where varargin is a lit of paired arguments. 
			%      arguments can be: 'analysis_type': 'stressbalance','thermal','masstransport','transient'
			%
			%   Examples:
			%      self=solve(self,'analysis_type','stressbalance');

			%recover options
			options=pairoptions(varargin{:});

			%add default options
			options=process_solve_options(options);

			%length of list
			nummodels=length(self.models);

			%name of queue: to make it unique, add a time stamp
			name=[self.name '-' datestr(now,1) '-' datestr(now,'HH-MM-SS') ];

			%name of cluster will be first name of list
			cluster=self.cluster;

			%Figure out parameters for this particular cluster
			cluster_rc_location=which('cluster.rc');
			[codepath,executionpath]=ClusterParameters(cluster,cluster_rc_location);

			%solve in batch mode: 
			for i=1:nummodels,

				%model
				mdex=self.models{i};

				%recover some fields
				mdex.analysis_type=options.analysis_type;

				mdex.name=[name '-' num2str(i) 'vs' num2str(nummodels)];
				mdex.time=self.time;
				mdex.queue=self.queue;
				mdex.cluster=self.cluster;
				if ~isnan(self.np),
					mdex.np=self.np;
				end

				%call solve in batch mode:
				if strcmpi(cluster,oshostname),
					mdex=solve(mdex,varargin{:});
				else
					mdex=solve(mdex,varargin{:},'batch','yes','directory',name);
				end

				%feed back
				self.models{i}=mdex;
			end

			%locally, we are done.
			if strcmpi(cluster,oshostname),
				return
			end

			%now, tar all the files and then erase them.
			eval(['!find -iname ''' name '-*'' > file_list.txt']);
			!tar zcvf ModelList.tar.gz --files-from file_list.txt
			!rm -rf *.bin *.queue file_list.txt

			%still have to build a launching script.
			BuildMultipleQueueingScript(cluster,name,executionpath,codepath);

			%launch jobs on remote cluster
			LaunchMultipleQueueJob(cluster,name,executionpath);

			%erase files: 
			delete([name '.queue']);
			delete('ModelList.tar.gz');

			%save name: 
			self.name=name;
		end % }}}
	end
end

function BuildMultipleQueueingScript(cluster,name,executionpath,codepath)% {{{
%BUILDMULTIPLEQUEUEINGSCRIPT - 
%
%   Usage:
%      BuildMultipleQueueingScript(executionpath,codepath)

disp('building queuing script');

%First try and figure out if there is a special script for this particular cluster
function_name=['BuildMultipleQueueingScript' cluster]

%some specific treatment of identical cluster, gemini, castor and pollux
if strcmpi(cluster,'castor') || strcmpi(cluster,'pollux'),
	function_name='BuildMultipleQueueingScriptgemini';
end

if exist(function_name,'file'),
	%Call this function:
	eval([function_name '(name,executionpath,codepath);']);
else
	%Call the generic BuildQueueingScript:
	BuildMultipleQueueingScriptGeneric(name,executionpath,codepath);
end
end % }}}
function BuildQueueingScriptgemini(name,executionpath,codepath)% {{{
%BUILDQUEUEINGSCRIPTGEMINI - ...
%
%   Usage:
%      BuildQueueingScriptgemini(md,executionpath,codepath)

scriptname=[name '.queue'];

fid=fopen(scriptname,'w');
if fid==-1,
	error(['BuildQueueingScriptgeminierror message: could not open ' scriptname ' file for ascii writing']);
end

fprintf(fid,'#!/bin/sh\n');
fprintf(fid,'cd %s\n',executionpath);
fprintf(fid,'mkdir %s\n',name);
fprintf(fid,'cd %s\n',name);
fprintf(fid,'mv ../ModelList.tar.gz ./\n');
fprintf(fid,'tar -zxvf ModelList.tar.gz\n');
fprintf(fid,'foreach i (%s-*vs*.queue)\n',name);
fprintf(fid,'qsub $i\n');
fprintf(fid,'end\n');
fclose(fid);
end% }}}
function LaunchMultipleQueueJob(cluster,name,executionpath)% {{{
%LAUNCHMULTIPLEQUEUEJOB - ...
%
%   Usage:
%      LaunchMultipleQueueJob(executionpath)

%First try and figure out if there is a special script for thie particular cluster
function_name=['LaunchMultipleQueueJob' cluster]

%some specific treatment of identical cluster, gemini, castor and pollux
if strcmpi(cluster,'castor') || strcmpi(cluster,'pollux'),
	function_name='LaunchMultipleQueueJobgemini';
end

if exist(function_name,'file'),
	%Call this function:
	eval([function_name '(cluster,name,executionpath);']);
else
	%Call the generic LaunchMultipleQueueJob:
	LaunchMultipleQueueJobGeneric(cluster,name,executionpath);
end
end% }}}
function md=LaunchMultipleQueueJobgemini(cluster,name,executionpath)% {{{
%LAUNCHMULTIPLEQUEUEJOBGEMINI - Launch multiple queuing script on Gemini cluster
%
%   Usage:
%      LaunchMultipleQueueJobgemini(cluster,name,executionpath)

%first, check we have the binary file and the queuing script
if ~exist([ name '.queue'],'file'),
	error('LaunchMultipleQueueJobgemini error message: queuing script issing, cannot go forward');
end

if ~exist('ModelList.tar.gz','file'),
	error('LaunchMultipleQueueJobgemini error message: inputs models file missing, cannot go forward');
end

%upload both files to cluster
disp('uploading input file,  queuing script and variables script');
system(['scp ModelList.tar.gz ' name '.queue '  cluster ':' executionpath]);

disp('launching solution sequence on remote cluster');
issmssh(cluster,login,['"cd ' executionpath ' && source ' name '.queue "']);
end% }}}
