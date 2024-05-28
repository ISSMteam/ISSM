%
%  function to write response levels
%
function []=rlev_write(fidi,dresp,params)

if isempty(dresp)
    return;
end

%  put responses into lists for writing

nresp=0;
respl={};
probl={};
rell ={};
grell={};

fnames=fieldnames(dresp);
for i=1:numel(fnames)
    nresp=nresp+numel(dresp.(fnames{i}));
    [respli,probli,relli,grelli]=prop_levels(dresp.(fnames{i}));
    respl=[respl respli];
    probl=[probl probli];
    rell =[rell  relli ];
    grell=[grell grelli];
end

%  write response levels

param_write(fidi,'\t  ','distribution',' ','\n',params);
if ~isempty(respl)
    rlevi_write(fidi,'response_levels',respl);
    param_write(fidi,'\t  ','compute',' ','\n',params);
end 
if ~isempty(probl)
    rlevi_write(fidi,'probability_levels',probl);
end
if ~isempty(rell)
    rlevi_write(fidi,'reliability_levels',rell);
end
if ~isempty(grell)
    rlevi_write(fidi,'gen_reliability_levels',grell);
end

end

%
%  function to each type of response level
%
function []=rlevi_write(fidi,ltype,levels)

fprintf(fidi,'\t  num_%s =',ltype);
for i=1:numel(levels)
    fprintf(fidi,' %d',length(levels{i}));
end
fprintf(fidi,'\n');

fprintf(fidi,'\t  %s =\n',ltype);

for i=1:numel(levels)
    if ~isempty(levels{i})
        vector_write(fidi,sprintf('\t    '),levels{i},8,76);
    end
end

end
