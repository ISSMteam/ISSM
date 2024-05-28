%
%  function to find a number of elements from a list of arrays.
%  
%  [asize]=array_numel(varargin)
%
%  see array_size to check the number and shape of elements, if
%  multiple indices will be used.
%
function [anum]=array_numel(varargin)

anum=1;

for iarg=1:nargin
    if ischar(varargin{iarg})
        inum=numel(cellstr(varargin{iarg}));
    else
        inum=numel(varargin{iarg});
    end

    if ~isequal(inum,1)
        if isequal(anum,1)
            anum=inum;
        else
            if ~isequal(inum,anum)
                if ~isempty(inputname(iarg))
                    error('Array ''%s'' has inconsistent number of elements.',inputname(iarg));
                else
                    error('Array %d has inconsistent number of elements.',iarg);
                end
            end
        end
    end
end

end
