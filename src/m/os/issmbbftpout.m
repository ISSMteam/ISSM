function issmbbftpout(host,path,login,port,numstreams,packages)
%BBFTPOUT put packages onto host, using bbftp. assuming unix system here.
%
%   usage: bbftpout(host,path,login,port,numstream,packages)
%
%

%get hostname
hostname=oshostname();

%if hostname and host are the same, do a simple copy
if strcmpi(host,hostname),
	for i=1:numel(packages),
		here=pwd;
		eval(['cd ' path])
		system(['rm -rf ' packages{i} ]);
		system(['ln -s ' here '/' packages{i} ' .']);
		eval(['cd ' here]);
	end
else 

	%build a string of the type: bbftp -s -u elarour -e 'setnbstream 8; cd /nobackupp10/elarour/Testing/Interactive3/; put Antarctica.tar.gz' pfe1.nas.nasa.gov
	command=['!bbftp -s -V -u ' login ' -e ''setnbstream 8; cd ' path '; '];
	for i=1:length(packages),
		command=[command 'put ' packages{i} ';'];
	end
	command=[command '''  pfe22.nas.nasa.gov'];

	eval(command);
end
