function expflip(domainname)
%EXPFLIP - flip orientation of all contours and domains in domainname exp file
%
%   Usage:
%      expflip(domainname)
%
%   Example:
%      expflip('MassFlux1.exp')

a=expread(domainname);

for i=1:length(a)
	a(i).x=flipud(a(i).x);
	a(i).y=flipud(a(i).y);
end

expwrite(a,domainname);
