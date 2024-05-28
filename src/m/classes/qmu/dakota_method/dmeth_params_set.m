%
%  set parameters of a dakota_method object.
%
%  [dm]=dmeth_params_set(dm,varargin)
%
function [dm]=dmeth_params_set(dm,varargin)

if ~isa(dm,'dakota_method')
    error('Object ''%s'' is a ''%s'' class object, not ''%s''.',...
        inputname(1),class(dm),'dakota_method');
end

%  loop through each parameter field in the input list

for i=1:2:length(varargin)
    if isfield(dm.params,varargin{i})
        dm.params.(varargin{i})=varargin{i+1};
    else
        warning('dmeth_params_set:unknown_param',...
            'No parameter ''%s'' for dakota_method ''%s''.',...
            varargin{i},dm.method);
    end
end

end
