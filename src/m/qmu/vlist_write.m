%
%  function to write variable list
%
function []=vlist_write(fidi,cstring,cstring2,dvar)

if isempty(dvar)
    return;
end

fnames=fieldnames(dvar);

%  put variables into lists for writing

nvar=0;
pinitpt=[];
plower =[];
pupper =[];
pmean  =[];
pstddev=[];
pinitst=[];
pstype =[];
pscale =[];
pabscissas =[];
ppairs_per_variable =[];
pcounts=[];
pdesc  =[];

for i=1:numel(fnames)
    nvar=nvar+numel(dvar.(fnames{i}));
    pinitpt=[pinitpt prop_initpt(dvar.(fnames{i}))];
    plower =[plower  prop_lower(dvar.(fnames{i})) ];
    pupper =[pupper  prop_upper(dvar.(fnames{i})) ];
    pmean  =[pmean   prop_mean(dvar.(fnames{i}))  ];
    pstddev=[pstddev prop_stddev(dvar.(fnames{i}))];
    pinitst=[pinitst prop_initst(dvar.(fnames{i}))];
    pstype =[pstype  prop_stype(dvar.(fnames{i})) ];
    pscale =[pscale  prop_scale(dvar.(fnames{i})) ];
    ppairs_per_variable =[ppairs_per_variable  prop_pairs_per_variable(dvar.(fnames{i})) ];
    pabscissas =[pabscissas  prop_abscissas(dvar.(fnames{i})) ];
    pcounts =[pcounts  prop_counts(dvar.(fnames{i})) ];
    pdesc  =[pdesc   prop_desc(dvar.(fnames{i}),fnames{i})];
end

%  write variables
%  (using Dakota 4.1 syntax for backward compatability)

disp(sprintf('  Writing %d %s variables.',nvar,cstring));

fprintf(fidi,'\t%s = %d\n',cstring,nvar);
if ~isempty(pinitpt)
    fprintf(fidi,'\t  %s_initial_point =\n',cstring2);
    vector_write(fidi,sprintf('\t    '),pinitpt,6,76);
end
if ~isempty(plower)
    fprintf(fidi,'\t  %s_lower_bounds =\n',cstring2);
    vector_write(fidi,sprintf('\t    '),plower ,6,76);
end
if ~isempty(pupper)
    fprintf(fidi,'\t  %s_upper_bounds =\n',cstring2);
    vector_write(fidi,sprintf('\t    '),pupper ,6,76);
end
if ~isempty(pmean)
    fprintf(fidi,'\t  %s_means =\n',cstring2);
    vector_write(fidi,sprintf('\t    '),pmean  ,6,76);
end
if ~isempty(pstddev)
    fprintf(fidi,'\t  %s_std_deviations =\n',cstring2);
    vector_write(fidi,sprintf('\t    '),pstddev,6,76);
end
if ~isempty(pinitst)
    fprintf(fidi,'\t  %s_initial_state =\n',cstring2);
    vector_write(fidi,sprintf('\t    '),pinitst,6,76);
end
if ~isempty(pstype)
    fprintf(fidi,'\t  %s_scale_types =\n',cstring2);
    vector_write(fidi,sprintf('\t    '),pstype ,6,76);
end
if ~isempty(pscale)
    fprintf(fidi,'\t  %s_scales =\n',cstring2);
    vector_write(fidi,sprintf('\t    '),pscale ,6,76);
end
if ~isempty(ppairs_per_variable)
    %fprintf(fidi,'\t  %s_pairs_per_variable =\n',cstring2);
    fprintf(fidi,'\t  pairs_per_variable =\n');
    vector_write(fidi,sprintf('\t    '),ppairs_per_variable ,6,76);
end
if ~isempty(pabscissas)
    %fprintf(fidi,'\t  %s_abscissas =\n',cstring2);
    fprintf(fidi,'\t  abscissas =\n');
    vector_write(fidi,sprintf('\t    '),pabscissas ,6,76);
end
if ~isempty(pcounts)
    %fprintf(fidi,'\t  %s_counts =\n',cstring2);
    fprintf(fidi,'\t  counts =\n');
    vector_write(fidi,sprintf('\t    '),pcounts ,6,76);
end
if ~isempty(pdesc)
    fprintf(fidi,'\t  descriptors =\n');
    vector_write(fidi,sprintf('\t    '),pdesc  ,6,76);
end

end
