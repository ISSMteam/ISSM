%
%  function to return the string size of an array
%
%  function [ssize]=string_size(a,varargin)
%
function [ssize]=string_size(a,varargin)

if ~nargin
    help string_size
    return
end

%  check for column or row vector (Matlab uses a minimum of two
%  dimensions, so this won't match Matlab standard output)

for iarg=1:nargin-1
    if strcmpi(varargin{iarg},'vector')
        if (ndims(a) == 2) && ((size(a,1) == 1) || (size(a,2) == 1))
            ssize =['(' num2str(numel(a)) ')'];
            return
        end
    end
end

%  do the general case

asize=size(a);

%  assemble string for output

ssize ='(';
for i=1:length(asize)-1;
    ssize =[ssize num2str(asize(i)) 'x'];
end
ssize =[ssize num2str(asize(end)) ')'];

end

