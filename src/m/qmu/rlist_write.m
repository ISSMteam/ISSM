%
%  function to write response list
%
function [rdesc]=rlist_write(fidi,cstring,cstring2,dresp,rdesc)

if isempty(dresp)
    return;
end

%  put responses into lists for writing
%  (and accumulate descriptors into list for subsequent writing)

nresp=0;
pstype =[];
pscale =[];
pweight=[];
plower =[];
pupper =[];
ptarget=[];

fnames=fieldnames(dresp);
for i=1:numel(fnames)
    nresp=nresp+numel(dresp.(fnames{i}));
    pstype =[pstype  prop_stype(dresp.(fnames{i})) ];
    pscale =[pscale  prop_scale(dresp.(fnames{i})) ];
    pweight=[pweight prop_weight(dresp.(fnames{i}))];
    plower =[plower  prop_lower(dresp.(fnames{i})) ];
    pupper =[pupper  prop_upper(dresp.(fnames{i})) ];
    ptarget=[ptarget prop_target(dresp.(fnames{i}))];
    rdesc  =[rdesc   prop_desc(dresp.(fnames{i}),fnames{i})];
end

%  write responses

disp(sprintf('  Writing %d %s responses.',nresp,cstring));

if strcmp(cstring,'calibration_terms')==1
	fprintf(fidi,'\t%s = %d\n',cstring,nresp);
	
else
	fprintf(fidi,'\tnum_%s = %d\n',cstring,nresp);
end

if ~isempty(pstype)
fprintf(fidi,'\t  %s_scale_types =\n',cstring2);
vector_write(fidi,sprintf('\t    '),pstype ,6,76);
end

if ~isempty(pscale)
    fprintf(fidi,'\t  %s_scales =\n',cstring2);
    vector_write(fidi,sprintf('\t    '),pscale ,6,76);
end
if ~isempty(pweight)
    switch cstring2
        case 'objective_function'
            fprintf(fidi,'\t  %s_weights =\n','multi_objective');
            vector_write(fidi,sprintf('\t    '),pweight,6,76);
        case 'least_squares_term'
            fprintf(fidi,'\t  %s_weights =\n','least_squares');
            vector_write(fidi,sprintf('\t    '),pweight,6,76);
    end
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

end
