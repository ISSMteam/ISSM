function [B E]=pkriging(x,y,observations,x_interp,y_interp,varargin)
%PKRIGING - parallel Kriging
%
%   Usage:
%      [B E]=pkriging(x,y,observations,x_interp,y_interp,varargin);

options=pairoptions(varargin{:});
cluster=getfieldvalue(options,'cluster',generic('np',1));
options=removefield(options,'cluster',0);
name   = ['krig' num2str(feature('GetPid'))];

%Prepare directory in execution
if strcmpi(cluster.name, oshostname())
	localexecdir = cluster.executionpath;
else
	localexecdir = [issmdir() '/execution/'];
end
if ~exist(localexecdir, 'dir')
	error(['Could not find directory ' issmdir() '/execution/']);
elseif numel(dir(localexecdir))>200
	warning([localexecdir ' has more than 200 subdirectories. Consider cleaning up your execution directory'])
end
root = [localexecdir '/' name];
if exist(root, 'dir')
	rmdir(root, 's');
end
mkdir(root);
basename = [root '/' name];

% =========================================   MARSHALL.m =================================================
disp(['marshalling file ' name '.bin']);
fid=fopen([basename '.bin'],'wb');
if(fid==-1) error(['marshall error message: could not open ' name '.bin file for binary writing']); end

%Write all data
WriteData(fid,'','name','md.x','data',x,'format','DoubleMat');
WriteData(fid,'','name','md.y','data',y,'format','DoubleMat');
WriteData(fid,'','name','md.data','data',observations,'format','DoubleMat');
WriteData(fid,'','name','md.x_interp','data',x_interp,'format','DoubleMat');
WriteData(fid,'','name','md.y_interp','data',y_interp,'format','DoubleMat');
options.marshall(fid);
WriteData(fid,'','name','md.EOF','data',true,'format','Boolean');
fclose(fid);

%Fake md as a place holder
md=model; md.cluster=cluster; md.settings.waitonlock=Inf; md.private.runtimename=name;md.miscellaneous.name=name;

%Launch job on remote cluster
BuildQueueScript(cluster, md, [basename '.queue'], 'kriging.exe');
UploadQueueJob(cluster,name,name,{[basename '.bin'] [basename '.queue']})
LaunchQueueJob(cluster,name,name,{[basename '.bin'] [basename '.queue']},'',0);

%Call waitonlock
waitonlock(md);

%Download
Download(cluster,name,{[name '.outbin']});
structure=parseresultsfromdisk(md,[name '.outbin'],0);
delete([name '.outlog']);
delete([name '.errlog']);
delete([name '.outbin']);
delete([name '.bin']);
delete([name '.tar.gz']);

%Process results
B=structure.predictions;
B=reshape(B,size(x_interp,2),size(x_interp,1))';
E=structure.error;
E=reshape(E,size(x_interp,2),size(x_interp,1))';
