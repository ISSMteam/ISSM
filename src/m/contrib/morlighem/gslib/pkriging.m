function [B E]=pkriging(x,y,observations,x_interp,y_interp,varargin)
%PKRIGING - parallel Kriging
%
%   Usage:
%      [B E]=pkriging(x,y,observations,x_interp,y_interp,varargin);

options=pairoptions(varargin{:});
cluster=getfieldvalue(options,'cluster',generic('np',1));
options=removefield(options,'cluster',0);
name   = ['krig' num2str(feature('GetPid'))];

if 1,
% =========================================   MARSHALL.m =================================================
disp(['marshalling file ' name '.bin']);
fid=fopen([name '.bin'],'wb');
if fid==-1,
	error(['marshall error message: could not open ' name '.bin file for binary writing']);
end

%Write all data
WriteData(fid,'','name','md.x','data',x,'format','DoubleMat');
WriteData(fid,'','name','md.y','data',y,'format','DoubleMat');
WriteData(fid,'','name','md.data','data',observations,'format','DoubleMat');
WriteData(fid,'','name','md.x_interp','data',x_interp,'format','DoubleMat');
WriteData(fid,'','name','md.y_interp','data',y_interp,'format','DoubleMat');

%Now, write number of options
options.marshall(fid);

%Last, write "md.EOF" to make sure that the binary file is not corrupt
WriteData(fid,'','name','md.EOF','data',true,'format','Boolean');

%Launch job on remote cluster
BuildKrigingQueueScript(cluster,name,'',1,0,0); %gather, valgrind, gprof
UploadQueueJob(cluster,name,name,{[name '.bin'] [name '.queue']})
LaunchQueueJob(cluster,name,name,{[name '.bin'] [name '.queue']},'',0);

%Call waitonlock
md=model; md.cluster=cluster; md.settings.waitonlock=Inf; md.private.runtimename=name;md.miscellaneous.name=name;
waitonlock(md);

%Download
end
Download(cluster,name,{[name '.outbin']});
structure=parseresultsfromdisk(md,[name '.outbin'],0);
delete([name '.outlog']);
delete([name '.errlog']);
delete([name '.outbin']);
delete([name '.bin']);
if ~ispc(),
	delete([name '.tar.gz']);
end

%Process results
B=structure.predictions;
B=reshape(B,size(x_interp,2),size(x_interp,1))';
E=structure.error;
E=reshape(E,size(x_interp,2),size(x_interp,1))';
