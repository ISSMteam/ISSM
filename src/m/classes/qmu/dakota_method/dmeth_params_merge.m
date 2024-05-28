%
%  merge a structure of parameters into a dakota_method object.
%
%  [dm]=dmeth_params_merge(dm,params)
%
function [dm]=dmeth_params_merge(dm,params)

if ~isa(dm,'dakota_method')
    error('Object ''%s'' is a ''%s'' class object, not ''%s''.',...
        inputname(1),class(dm),'dakota_method');
end

%  loop through each parameter field in the structure

fnames=fieldnames(params);

for i=1:numel(fnames)
    if isfield(dm.params,fnames{i})
        dm.params.(fnames{i})=params.(fnames{i});
    else
        warning('dmeth_params_merge:unknown_param',...
            'No parameter ''%s'' for dakota_method ''%s''.',...
            fnames{i},dm.method);
    end
end

end
