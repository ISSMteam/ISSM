%
%  function to convert an item to a string
%
%  function [svec]=item2str(a)
%
function [svec]=item2str(a)

if     islogical(a)
    if a
        svec='true';
    else
        svec='false';
    end
elseif ischar(a)
    svec=['''' a ''''];
elseif isnumeric(a)
    svec=num2str(a);
else
    if ~isempty(inputname(1))
        warning('item2str:item_unrecog',...
            'Item ''%s'' is of unrecognized type ''%s''.',...
            inputname(1),class(a));
    else
        warning('item2str:item_unrecog',...
            'Item %d is of unrecognized type ''%s''.',...
            1,class(a));
    end
    return
end

end
