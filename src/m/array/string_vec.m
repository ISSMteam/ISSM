function [svec]=string_vec(a)
%STRING_VEC - return the string of a vector
%
%   Usage:
%      svec=string_vec(a);

if ~nargin
    help string_vec
    return
end

if (numel(a) == 0)
    svec='[]';
    return
end

%  assemble string for output

svec ='[';
for i=1:numel(a)-1;
    svec=[svec item2str(a(i)) ' '];
end
svec=[svec item2str(a(end)) ']'];

end
