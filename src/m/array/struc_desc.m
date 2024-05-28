%
%  function to find the structures with the specified descriptors
%  
%  [sarrayo]=struc_desc(sarray,varargin)
%
function [sarrayo]=struc_desc(sarray,varargin)

if ~isfield(sarray,'descriptor')
    if ~isempty(inputname(1))
        error('Field ''descriptor'' not found in array ''%s''.',inputname(1));
    else
        error('Field ''descriptor'' not found in array %d.',1);
    end
end

sarrayo=struct([]);

for iarg=1:nargin-1
    if     iscell(varargin{iarg})
        desc=        varargin{iarg};
    elseif ischar(varargin{iarg})
        desc=cellstr(varargin{iarg});
    end

    for i=1:length(desc)
        sarrayoi=struc_desci(sarray,desc{i});
        if ~isempty(sarrayoi)
            if isempty(sarrayo)
                sarrayo       =sarrayoi;
            else
                sarrayo(end+1)=sarrayoi;
            end
        end
    end
end

%  if nothing found, return whole array

if isempty(sarrayo)
    sarrayo=sarray;
end

end

%
%  function to find the structure with the specified descriptor
%  
function [sarrayo]=struc_desci(sarray,str)

sarrayo=struct([]);

for i=1:numel(sarray)
    if strcmp(sarray(i).descriptor,str)
        sarrayo=sarray(i);
        return
    end
end

warning(['String ''' str ''' not found in array ''' inputname(1) '''.']);

end
