function [svec]=any2str(a,alim)
%ANY2STR - convert anything to a string
%
%   Usage:
%      svec=any2str(a,alim);

if ~exist('alim','var') || (numel(a) <= alim)
    if iscell(a)
        svec=string_cell(a);
    else
        if (numel(a) > 1) && ~ischar(a)
            svec=string_vec(a);
        else
            svec=item2str(a);
        end
    end
else
	svec=[string_size(a) ' ''' class(a) ''''];
end

end
