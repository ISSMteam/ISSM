function issmssh(host,login,port,command)
%ISSMSSH - wrapper for OS independent ssh command.
%
%   usage: 
%      issmssh(host,command)

%first get hostname 
hostname=oshostname();

%if same as host, just run the command. 
if strcmpi(host,hostname),
	system(command);
else
	if ispc & ~ismingw,
		%use the putty project plink.exe: it should be in the path.

		%get ISSM_DIR variable
		[status,ISSM_DIR]=system('echo [%ISSM_DIR_WIN%]');
		if status, 
			error('issmssh error message: could not find ISSM_DIR_WIN environment variable');
		end
		ISSM_DIR=ISSM_DIR(2:end-2);

		username=input('Username: (quoted string) ');
		key=input('Key: (quoted string) ');

		system([ISSM_DIR '/externalpackages/ssh/plink.exe -ssh -l "' username '" -pw "' key '" ' host ' "' command '"']);

	else
		%just use standard unix ssh
		if port,
			eval(['!ssh -l ' login ' -p ' num2str(port) ' localhost "' command '"']);
		else
			eval(['!ssh -l ' login ' ' host ' "' command '"']);
		end
	end
end
