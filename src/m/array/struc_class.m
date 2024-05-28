%
%  function to find the structural fields of a specified class
%  
%  [sclasso]=struc_class(sclass,cstr)
%
function [sclasso]=struc_class(sclass,cstr)

%  collect only the objects of the appropriate class

if     isa(sclass,cstr)
    if ~isempty(inputname(1))
        sclasso.(inputname(1))=sclass;
    else
        sclasso.(cstr)        =sclass;
    end

elseif isstruct(sclass)
    fnames=fieldnames(sclass);
    for i=1:numel(fnames)
        if isa(sclass.(fnames{i}),cstr)
            sclasso.(fnames{i})=sclass.(fnames{i});
        end
    end
end

if ~exist('sclasso','var')
    sclasso=struct([]);
end

end
