%
%  function to return the string of an vector
%
%  function [svec]=string_vec(a)
%
function [svec]=string_vec(a)

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
