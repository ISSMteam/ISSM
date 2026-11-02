function [anum]=array_numel(varargin)
%ARRAY_NUMEL - find the number of elements from a list of arrays
%
%   See array_size to check the number and shape of elements, if
%   multiple indices will be used.
%
%   Usage:
%      anum=array_numel(a1,a2,...);

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
