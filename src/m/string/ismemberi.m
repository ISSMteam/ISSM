function flag=ismemberi(string,list),
%ISMEMBERI - return 1 if a string belongs to a list (case insensitive)
%
%   same function as Matlab's ismember except that it
%   is case insensitive
%
%   Usage:
%      flag=ismemberi(string,list);
%
%   Example:
%      flag=ismemberi('test','{'test1','test2','test3'});

if ~iscell(list)
	error('ismemberi error message: the list of string must be a cell!')
end

%initialize output
flag=0;

%go through the list
for i=1:length(list),
	if strcmpi(string,list{i}),
		flag=i;
		return
	end
end
