function issmscpout(host,path,login,port,packages,varargin)
%ISSMSCPOUT send packages to a host, using scp on unix, and pscp on windows
%
%   usage: issmscpout(host,path,login,port,packages)
%
%

%get hostname
hostname=oshostname();

%are we disallowing symbolic links? 
if nargin==6,
	no_symlinks=1;
else
	no_symlinks=0;
end

%if hostname and host are the same, do a simple copy

if strcmpi(host,hostname),
	for i=1:numel(packages),
		here=pwd;
		eval(['cd ' path])
		system(['rm -rf ' packages{i} ]);
		if no_symlinks,
			system(['cp ' here '/' packages{i} ' .']);
		else
			system(['ln -s ' here '/' packages{i} ' .']);
		end
		eval(['cd ' here]);
	end
else 
	if ispc & ~ismingw,
		%use the putty project pscp.exe: it should be in the path.

		%get ISSM_DIR variable
		[status,ISSM_DIR]=system('echo [%ISSM_DIR_WIN%]');
		if status, 
			error('issmscpout error message: could not find ISSM_DIR_WIN environment variable');
		end
		ISSM_DIR=ISSM_DIR(2:end-2);

		username=input('Username: (quoted string) ');
		key=input('Key: (quoted string) ');

		for i=1:numel(packages),
			[status,result]=system([ISSM_DIR '/externalpackages/ssh/pscp.exe -l "' username '" -pw "' key '" ' packages{i} ' ' host ':' path]);
			if status, 
				error('issmscpout error message: could not call putty pscp');
			end
		end

	else
		%just use standard unix scp
		%create string of packages being sent
		string='';
		for i=1:numel(packages),
			string=[string ' ' packages{i}];
		end
		string=[string ' '];

		if port,
			eval(['!scp -P ' num2str(port) ' ' string ' ' login '@localhost:' path]);
		else
			eval(['!scp ' string ' ' login '@' host ':' path]);
		end
	end
end
