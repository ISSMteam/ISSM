%
%  function to write linear constraint list
%
function []=lclist_write(fidi,cstring,cstring2,dvar)

if isempty(dvar)
    return;
end

%  put linear constraints into lists for writing

nvar=0;
pmatrix=[];
plower =[];
pupper =[];
ptarget=[];
pstype =[];
pscale =[];

fnames=fieldnames(dvar);
for i=1:numel(fnames)
    nvar=nvar+numel(dvar.(fnames{i}));
    pmatrix=[pmatrix prop_matrix(dvar.(fnames{i}))];
    plower =[plower  prop_lower(dvar.(fnames{i})) ];
    pupper =[pupper  prop_upper(dvar.(fnames{i})) ];
    ptarget=[ptarget prop_target(dvar.(fnames{i}))];
    pstype =[pstype  prop_stype(dvar.(fnames{i})) ];
    pscale =[pscale  prop_scale(dvar.(fnames{i})) ];
end

%  write linear constraints

disp(sprintf('  Writing %d %s linear constraints.',...
    nvar,cstring));

if ~isempty(pmatrix)
    fprintf(fidi,'\t  %s_matrix =\n',cstring2);
    vector_write(fidi,sprintf('\t    '),pmatrix,6,76);
end
if ~isempty(plower)
    fprintf(fidi,'\t  %s_lower_bounds =\n',cstring2);
    vector_write(fidi,sprintf('\t    '),plower ,6,76);
end
if ~isempty(pupper)
    fprintf(fidi,'\t  %s_upper_bounds =\n',cstring2);
    vector_write(fidi,sprintf('\t    '),pupper ,6,76);
end
if ~isempty(ptarget)
    fprintf(fidi,'\t  %s_targets =\n',cstring2);
    vector_write(fidi,sprintf('\t    '),ptarget,6,76);
end
if ~isempty(pstype)
    fprintf(fidi,'\t  %s_scale_types =\n',cstring2);
    vector_write(fidi,sprintf('\t    '),pstype ,6,76);
end
if ~isempty(pscale)
    fprintf(fidi,'\t  %s_scales =\n',cstring2);
    vector_write(fidi,sprintf('\t    '),pscale ,6,76);
end

end
