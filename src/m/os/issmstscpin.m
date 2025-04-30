function issmstscpin(host, login,path, packages)
%ISSMSCPIN get packages from host, using scp on unix, and pscp on windows
%
%   usage: issmstscpin(host,long, path, packages)
%

%get initial warning mode
state=warning('query', 'all');
%remove warnings in case the files do not exist
warning off
for i=1:numel(packages),
	delete(packages{i});
end
%back to initial warning state
warning(state);

%use starcluster to scp 
%string to copy multiple files using scp: 
if numel(packages)==1,
	string=packages{1};
else
	string='';
	for i=1:numel(packages)-1,
		string=[string ' ' path packages{i} ' '];
	end
	string=[string path packages{end}];
end

system([starcluster() ' get ' host ' -u ' login ' ' string  ' ./']);

%check starcluster get  worked
for i=1:numel(packages),
	if ~exist(['./' packages{i}]),
		error('issmstscpin error message: could not call scp on *nix system');
	end
end
