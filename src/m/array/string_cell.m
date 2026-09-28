function [svec]=string_cell(a)
%STRING_CELL - return the string of a cell array
%
%   Usage:
%      svec=string_cell(a);

if ~nargin
    help string_cell
    return
end

if (numel(a) == 0)
    svec='{}';
    return
end

%  assemble string for output

svec ='{';
for i=1:numel(a)-1;
    svec=[svec item2str(a{i}) ' '];
end
svec=[svec item2str(a{end}) '}'];

end
