function list=listfiles()
%LISTFILES list files inside a directory
%        this is very OS dependent.
%
%   usage: list=listfiles;
%
%
%   see also LS DIR

%use dir, as it seems to act OS independent

first_list=dir;
list={};

for i=1:numel(first_list),
	if (  ~strcmpi(first_list(i).name,'.') &...
			~strcmpi(first_list(i).name,'..') &...
			~strcmpi(first_list(i).name,'NightlyRun') &...
			~strcmpi(first_list(i).name,'.svn')),
		list{end+1}=first_list(i).name;
	end
end
