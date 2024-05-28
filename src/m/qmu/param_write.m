%
%  function to write a parameter
%
function []=param_write(fidi,sbeg,pname,smid,send,params)

if ~isfield(params,pname)
    warning('param_write:param_not_found',...
        'Parameter ''%s'' not found in structure.',pname);
    return
end

if islogical(params.(pname)) && ~params.(pname)
    return
end

if     islogical(params.(pname))
    fprintf(fidi,[sbeg '%s' send],pname);
elseif ischar   (params.(pname))
    fprintf(fidi,[sbeg '%s' smid '%s' send],pname,params.(pname));
elseif isnumeric(params.(pname))
    fprintf(fidi,[sbeg '%s' smid '%g' send],pname,params.(pname));
end

end
