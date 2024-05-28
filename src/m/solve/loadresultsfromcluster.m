function md=loadresultsfromcluster(md,varargin)
%LOADRESULTSFROMCLUSTER - load results of solution sequence from cluster
%
%   Usage:
%      md=loadresultsfromcluster(md,varargin);
%
%      Options include: 'runtimename', 'nolog'
%
%   Example:
%     md=loadresultsfromcluster(md,'runtimename','test101-06-15-2021-13-24-18-4883');

%process options: 
options=pairoptions(varargin{:});
nolog=getfieldvalue(options,'nolog',0);
md.private.runtimename = getfieldvalue(options,'runtimename',md.private.runtimename);

%retrieve cluster, to be able to call its methods
cluster=md.cluster;

%Download outputs from the cluster
if ~nolog,
	filelist={[md.miscellaneous.name '.outlog'],[md.miscellaneous.name '.errlog']};
else
	filelist={};
end
if md.qmu.isdakota,
	filelist{end+1}=[md.miscellaneous.name '.qmu.err'];
	filelist{end+1}=[md.miscellaneous.name '.qmu.out'];
	if isfield(md.qmu.params,'tabular_graphics_data'),
		if md.qmu.params.tabular_graphics_data==true,
			filelist{end+1}='dakota_tabular.dat';
		end
	end
	if md.qmu.output & strcmpi(md.qmu.statistics.method(1).name,'None'),
		if strcmpi(md.qmu.method.method,'nond_sampling'),
			for i=1:md.qmu.method.params.samples
				filelist{end+1}=[md.miscellaneous.name '.outbin.' num2str(i)];
			end
		end
	end
	if ~strcmpi(md.qmu.statistics.method(1).name,'None'),
		filelist{end+1}=[md.miscellaneous.name '.stats'];
	end
else
	filelist{end+1}=[md.miscellaneous.name '.outbin'];
end
Download(cluster,md.private.runtimename,filelist);

%If we are here, no errors in the solution sequence, call loadresultsfromdisk.
md=loadresultsfromdisk(md,[md.miscellaneous.name '.outbin']);

%erase the log and output files
for i=1:numel(filelist)
	filename = filelist{i};
	if exist(filename)
		delete(filename)
	end
end
if exist([md.private.runtimename '.tar.gz']) & ~ispc(),
	delete([md.private.runtimename '.tar.gz']);
end

%erase input file if run was carried out on same platform.
hostname=oshostname();
if strcmpi(hostname,cluster.name),
	delete([md.miscellaneous.name '.bin']);
	delete([md.miscellaneous.name '.toolkits']);
	if md.qmu.isdakota
		delete([md.miscellaneous.name '.qmu.in']);
	end
	if ~ispc(),
		delete([md.miscellaneous.name '.queue']);
	else
		delete([md.miscellaneous.name '.bat']);
	end
end
