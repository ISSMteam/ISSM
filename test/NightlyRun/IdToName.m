function name=IdToName(id);
%IDTONAME- return name of test
%  
%   if id=0, the full list of test is returned
%
%   Usage:
%      name=IdToName();
%      name=IdToName(id);

if nargin==0 | id==0,
	flist=dir;%use dir, as it seems to act OS independent
	list_ids=[];
	for i=1:numel(flist),
		if ( strncmp(flist(i).name,'test',4) &...                         %File name must start with 'test'
				strncmp(fliplr(flist(i).name),fliplr('.m'),2)&...           %File name must end by '.m'
				~strcmp(flist(i).name,'test.m'))                            %File name must be different than 'test.m'
			id=str2num(flist(i).name(5:end-2));
			if isempty(id),
				disp(['WARNING: ignore file ' flist(i).name ]);
			else
				list_ids(end+1)=eval(flist(i).name(5:end-2));                  %Keep test id only (skip 'test' and '.m')
			end
		end
	end
	list_ids=sort(list_ids);
	for i=list_ids,
		name=IdToName(i);
		disp(['test ' num2str(i,'%5i\n') ' : ' name]);
	end
	return;
end

filename = ['test' num2str(id) '.m'];

if ~exist(filename,'file')
	error(['file ' filename ' does not exist']);
end

string='%TestName:';
fid=fopen(filename,'r');
A=fscanf(fid,'%s',3);
if ~strncmp(A,string,numel(string)) | numel(A)<numel(string)+2,
	error(['Test file ' filename ' does to start with a test name']);
end
name = A(numel(string)+1:end);
fclose(fid);
