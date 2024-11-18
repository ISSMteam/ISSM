function issmscpin(host, login,port,path, packages)
%ISSMSCPIN get packages from host, using scp on unix, and pscp on windows
%
%   usage: issmscpin(host,packages,path)
%
%

%first get hostname
hostname=oshostname();

%first be sure packages are not in the current directory, this could conflict with pscp on windows. 
for i=1:numel(packages),
	if exist(packages{i},'file')
		delete(packages{i});
	end
end

%if hostname and host are the same, do a simple copy
if strcmpi(hostname,host)
	for i=1:numel(packages)
		system(['cp ' path '/' packages{i} ' .']);
	end
else
	%just use standard unix scp string to copy multiple files using scp: 
	if numel(packages)==1,
		string=packages{1};
	else
		string='\{';
		for i=1:numel(packages)-1,
			string=[string packages{i} ','];
		end
		string=[string packages{end} '\}'];
	end

	if port,
		eval(['!scp -P ' num2str(port) ' ' login '@localhost:' path '/' string ' ./']);
	else
		eval(['!scp ' login '@' host ':' path '/' string ' ./']);
	end

	%check scp worked
	for i=1:numel(packages),
		if ~exist(['./' packages{i}]),
			warning(['issmscpin error message: could not scp ' packages{i}]);
		end
	end
end
