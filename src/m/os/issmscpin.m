function issmscpin(host, login, port, path, packages, bracketstyle)
%ISSMSCPIN get packages from host
%
%   usage:
%      issmscpin(host,packages,path,bracketstyle)
%
%      bracketstyle:  1 - \{\}    (escaped; default)
%                     2 - {}      (not escaped)

%does machine require escaped brackets?
if nargin==5
	bracketstyle = 1;
end

%get hostname
hostname=oshostname();

%if hostname and host are the same, do a simple copy
if strcmpi(hostname,host)
	for i=1:numel(packages)
		system(['cp ' path '/' packages{i} ' .']);
	end
else
	if numel(packages)==1
		fileliststr=packages{1};
	else
		fileliststr='\{';
		for i=1:numel(packages)-1
			fileliststr=[fileliststr packages{i} ','];
		end
		fileliststr=[fileliststr packages{end} '\}'];

		%remove \ if bracketstyle is 2
		if bracketstyle==2
			fileliststr = [fileliststr(2:end-2) '}'];
		end
	end

	if port
		[status]=system(['scp -P ' num2str(port) ' ' login '@localhost:' path '/' fileliststr ' ./']);
		if status ~= 0
			%List expansion is a bashism. Try again with '-OT'.
			[status,cmdout]=system(['scp -OT -P ' num2str(port) ' ' login '@localhost:' path '/' fileliststr ' ./']);
		end
	else
		[status]=system(['scp ' login '@' host ':' path '/' fileliststr ' ./']);
		if status ~= 0
			%List expansion is a bashism. Try again with '-OT'.
			[status,cmdout]=system(['scp -OT ' login '@' host ':' path '/' fileliststr ' ./']);
		end
	end

	%check scp worked
	if status ~= 0
		error(['issmscpin error message: ' cmdout])
	end
	for i=1:numel(packages),
		if ~exist(['./' packages{i}]),
			warning(['issmscpin error message: could not scp ' packages{i}]);
		end
	end
end
