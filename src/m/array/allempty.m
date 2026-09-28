function [cout]=allempty(cin)
%ALLEMPTY - return an empty cell array if all array elements are empty
%
%   Usage:
%      cout=allempty(cin);

if ~nargin
    help allempty
    return
end

for j=1:numel(cin)
    if ~isempty(cin{j})
        cout=cin;
        return
    end
end
cout={};

end

