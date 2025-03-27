function issmscpin(host, login,port,path, packages)
%ISSMSCPIN get packages from host
%
%   usage: issmscpin(host,packages,path)
%

%first get hostname
hostname=oshostname();

%if hostname and host are the same, do a simple copy
if strcmpi(hostname,host)
	for i=1:numel(packages)
		system(['cp ' path '/' packages{i} ' .']);
	end
else
	if numel(packages)==1,
		fileliststr=packages{1};
	else
		fileliststr='\{';
		for i=1:numel(packages)-1,
			fileliststr=[fileliststr packages{i} ','];
		end
		fileliststr=[fileliststr packages{end} '\}'];
	end

	if port,
		[status,cmdout]=system(['scp -P ' num2str(port) ' ' login '@localhost:' path '/' fileliststr ' ./']);
		if status ~= 0,
			%List expansion is a bash'ism. Try again with -OT.
			[status,cmdout]=system(['scp -OT -P ' num2str(port) ' ' login '@localhost:' path '/' fileliststr ' ./']);
		end
	else
		[status,cmdout]=system(['scp ' login '@' host ':' path '/' fileliststr ' ./']);
		if status ~= 0,
			%List expansion is a bash'ism. Try again with -OT.
			[status,cmdout]=system(['scp -OT ' login '@' host ':' path '/' fileliststr ' ./']);
		end
	end

	%check scp worked
	if status ~= 0,
		error(['issmscpin error message: ' cmdout])
	end
	for i=1:numel(packages),
		if ~exist(['./' packages{i}]),
			warning(['issmscpin error message: could not scp ' packages{i}]);
		end
	end
end
