function issmstscpout(host,path,login,packages)
%ISSMSTSCPOUT send packages to a host, using starcluster put on unix
%
%   usage: issmstscpout(host,path,login,packages)
%
%

%create string of packages being sent
string='';
for i=1:numel(packages),
	string=[string ' ' packages{i}];
end
string=[string ' '];

system([ starcluster() ' put ' host ' -u ' login ' ' string ' ' path]);
