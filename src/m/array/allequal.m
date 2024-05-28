%  function to return an empty array if all array elements are
%  equal to the given value, which may also be empty but not nan.
%
%  (note that by definition, nan is not equal to nan.  this could
%  be changed by using isequalwithequalnans.)
%
%  function [aout]=allequal(ain,aval)
%
function [aout]=allequal(ain,aval)

if ~nargin
    help allequal
    return
end

aout=ain;

if     islogical(ain) && islogical(aval)
    for i=1:numel(ain)
        if ~isequal(ain(i),aval)
            return
        end
    end
    aout=logical([]);

elseif isnumeric(ain) && isnumeric(aval)
    for i=1:numel(ain)
        if ~isequal(ain(i),aval)
            return
        end
    end
    aout=[];

elseif ischar(ain) && ischar(aval)
    for i=1:size(ain,1)
        if ~strcmp(ain(i,:),aval)
            return
        end
    end
    aout='';

elseif iscell(ain)
    if     islogical(aval)
        for i=1:numel(ain)
            if ~islogical(ain{i}) || ~isequal(ain{i},aval)
                return
            end
        end
        aout={};

    elseif isnumeric(aval)
        for i=1:numel(ain)
            if ~isnumeric(ain{i}) || ~isequal(ain{i},aval)
                return
            end
        end
        aout={};

    elseif ischar(aval)
        for i=1:size(ain,1)
            if ~ischar(ain{i}) || ~strcmp(ain{i},aval)
                return
            end
        end
        aout={};
    end
end

end
