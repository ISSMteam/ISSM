function issmscpout(host, path, login, port, packages, no_symlinks, bracketstyle)
%ISSMSCPOUT send files to host
%
%   usage:
%      issmscpout(host,path,login,port,packages,no_symlinks,bracketstyle)
%
%      bracketstyle:  1 - \{\}    (escaped; default)
%                     2 - {}      (not escaped)

%get hostname
hostname=oshostname();

%disable symbolic links?
if nargin<6
	no_symlinks=0;
end
%does machine require escaped brackets?
if nargin<7
	bracketstyle = 1;
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

%General case: this is not a local machine
else
	if numel(packages)==1
		fileliststr=packages{1};
	else
		fileliststr='\{';
		for i=1:numel(packages)-1,
			fileliststr=[fileliststr packages{i} ','];
		end
		fileliststr=[fileliststr packages{end} '\}'];

		%remove backslashes if bracketstyle is 2
		if bracketstyle==2
			fileliststr = [fileliststr(2:end-2) fileliststr(end)];
		end
	end
	if port
		disp(['scp -P ' num2str(port) ' ' fileliststr ' ' login '@localhost:' path])
		[status]=system(['scp -P ' num2str(port) ' ' fileliststr ' ' login '@localhost:' path]);
		if status~=0
			%List expansion is a bashism. Try again with '-OT'.
			[status,cmdout]=system(['scp -OT -P ' num2str(port) ' ' fileliststr ' ' login '@localhost:' path]);
		end
	else
		[status]=system(['scp ' fileliststr ' ' login '@' host ':' path]);
		if status~=0
			%List expansion is a bashism. Try again with '-OT'.
			[status,cmdout]=system(['scp -OT ' fileliststr ' ' login '@' host ':' path]);
		end
	end

	%check scp worked
	if status~=0
		error(['issmscpin error message: ' cmdout])
	end
end
