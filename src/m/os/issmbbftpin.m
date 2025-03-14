function issmbbftpin(host, login,port,numstreams,path, packages)
%BBFTPIN get packages from host, using bbftp. assuming unix system here.
%
%   usage: scpin(host,packages,path)
%

%first get hostname
hostname=oshostname();

%first be sure packages are not in the current directory, this could conflict with pscp on windows. 
for i=1:numel(packages),
	if exist(packages{i},'file'),
		delete(packages{i});
	end
end

%if hostname and host are the same, do a simple copy
if strcmpi(hostname,host),

	for i=1:numel(packages),
		success=copyfile([path '/' packages{i}]); %keep going, even if success=0
	end

else

	%build a string of the type: bbftp -s -u elarour -e 'setnbstream 8; cd /nobackupp10/elarour/Testing/Interactive3/; get Antarctica.outbin' pfe1.nas.nasa.gov
	command=['!bbftp -s -V -u ' login ' -e ''setnbstream 8; cd ' path '; '];
	for i=1:numel(packages),
		command=[command 'get ' packages{i} ';'];
	end
	command=[command '''  pfe22.nas.nasa.gov'];

	eval(command);

	%check bbftp worked
	for i=1:numel(packages),
		if ~exist(['./' packages{i}],'file'),
			error('scpin error message: could not call scp on *nix system');
		end
	end
end
