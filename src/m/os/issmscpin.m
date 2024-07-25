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
	if ispc & ~ismingw,
		%use the putty project pscp.exe: it should be in the path.

		%get ISSM_DIR variable
		[status,ISSM_DIR]=system('echo [%ISSM_DIR_WIN%]');
		if status, 
			error('issmscpin error message: could not find ISSM_DIR_WIN environment variable');
		end
		ISSM_DIR=ISSM_DIR(2:end-2);

		username=input('Username: (quoted string) ');
		key=input('Key: (quoted string) ');

		for i=1:numel(packages),
			[status,result]=system([ISSM_DIR '/externalpackages/ssh/pscp.exe -l "' username '" -pw "' key '" ' host ':' path '/' packages{i} ' ./']);
			if status, 
				error('issmscpin error message: could not call putty pscp');
			end
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
end
