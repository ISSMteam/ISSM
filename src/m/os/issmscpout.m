function issmscpout(host,path,login,port,packages,varargin)
%ISSMSCPOUT send files to host
%
%   usage: issmscpout(host,path,login,port,packages)
%

%get hostname
hostname=oshostname();

%are we disallowing symbolic links? 
if nargin==6
	no_symlinks=1;
else
	no_symlinks=0;
end

%if hostname and host are the same, do a simple copy or symlinks
if strcmpi(host,hostname)

	%Process both paths and add \ if there are any white spaces
	here = replace(pwd(), ' ', '\ ');

	for i=1:numel(packages)
		system(['rm -rf ' path '/' packages{i} ]);
		if no_symlinks
			system(['cp ' packages{i} ' ' path]);
		else
			system(['ln -s ' here '/' packages{i} ' ' path]);
		end
	end

%General case, this is not a local machine
else
	if numel(packages)==1
		fileliststr=packages{1};
	else
		fileliststr='\{';
		for i=1:numel(packages)-1,
			fileliststr=[fileliststr packages{i} ','];
		end
		fileliststr=[fileliststr packages{end} '\}'];
	end
	if port
		disp(['scp -P ' num2str(port) ' ' fileliststr ' ' login '@localhost:' path])
		[status,cmdout]=system(['scp -P ' num2str(port) ' ' fileliststr ' ' login '@localhost:' path]);
		if status~=0
			%List expansion is a bash'ism. Try again with -OT.
			[status,cmdout]=system(['scp -OT -P ' num2str(port) ' ' fileliststr ' ' login '@localhost:' path]);
		end
	else
		[status,cmdout]=system(['scp ' fileliststr ' ' login '@' host ':' path]);
		if status~=0
			%List expansion is a bash'ism. Try again with -OT.
			[status,cmdout]=system(['scp -OT ' fileliststr ' ' login '@' host ':' path]);
		end
	end

	%check scp worked
	if status~=0
		error(['issmscpin error message: ' cmdout])
	end
end
