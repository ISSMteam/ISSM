%
%  function to find an array size from a list of arrays.
%  
%  [asize]=array_size(varargin)
%
%  see array_numel to check only the number of elements, if
%  single indices will be used.
%
function [asize]=array_size(varargin)

asize=[1 1];

for iarg=1:nargin
    if ischar(varargin{iarg})
        isize=size(cellstr(varargin{iarg}));
    else
        isize=size(varargin{iarg});
    end

    if ~isequal(isize,[1 1])
        if isequal(asize,[1 1])
            asize=isize;
        else
            if ~isequal(isize,asize)
                if ~isempty(inputname(iarg))
                    error('Array ''%s'' has inconsistent size.',inputname(iarg));
                else
                    error('Array %d has inconsistent size.',iarg);
                end
            end
        end
    end
end

end
