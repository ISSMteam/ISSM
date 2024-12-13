function issmscpin(host, login,port,path, packages)
%ISSMSCPIN get packages from host, using scp on unix, and pscp on windows
%
%   usage: issmscpin(host,packages,path)
%
% NOTE: If users again have issues with file list (i.e.
%
%   {<FILE1>,<FILE2>,...,<FILEN>}
%
% ), note that this a bash'ism and default shell should be checked. View file 
% history for potential fix (i.e. some combination of -O and -T options).
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
		eval(['!scp -P ' num2str(port) ' ' login '@localhost:' path '/' fileliststr ' ./']);
	else
		eval(['!scp ' login '@' host ':' path '/' fileliststr ' ./']);
	end

	%check scp worked
	for i=1:numel(packages),
		if ~exist(['./' packages{i}]),
			warning(['issmscpin error message: could not scp ' packages{i}]);
		end
	end
end
