%
%  read a Dakota .in input file and parse it.
%
%  [method,dvar,dresp]=dakota_in_parse(filei)
%
%  where the required input is:
%    filei         (character, name of .in file)
%
%  the required output is:
%    method        (character, dakota method name)
%    dvar          (structure array, variables)
%    dresp         (structure array, responses)
%
%  the filei will be prompted if empty.  the fields of dvar and
%  dresp are particular to the data contained within the file.
%
%  this function reads a dakota .in input file and parses it
%  into the matlab workspace.  it operates in a content-driven
%  fashion, where it parses whatever input data it encounters
%  in the file, rather than searching for data based on the
%  particular method.  (this makes it independent of method.)
%
%  as of now, parameters are generally not parsed.  also, the
%  variable and response classes are not used for output.
%
%  this data would typically be used for modifying and submitting
%  a subsequent dakota run.  it could also be used with output
%  data for post-processing or annotation purposes.
%
%  "Copyright 2009, by the California Institute of Technology.
%  ALL RIGHTS RESERVED. United States Government Sponsorship
%  acknowledged. Any commercial use must be negotiated with
%  the Office of Technology Transfer at the California Institute
%  of Technology.  (NTR 47078)
%
%  This software may be subject to U.S. export control laws.
%  By accepting this  software, the user agrees to comply with
%  all applicable U.S. export laws and regulations. User has the
%  responsibility to obtain export licenses, or other export
%  authority as may be required before exporting such information
%  to foreign countries or providing access to foreign persons."
%
function [method,dvar,dresp]=dakota_in_parse(filei)

if ~nargin
    help dakota_in_parse
    return
end

if ~exist('filei' ,'var') || isempty(filei)
    filei=input('Input file?  ','s');
end
fidi=fopen(sprintf('%s',filei),'r');
if (fidi < 0)
    error('%s could not be opened.',filei);
end

%%  loop through the file to find the Dakota method

method=[];
fseek(fidi,0,'bof');
[fline]=findline(fidi,'method');
if ~ischar(fline)
    return
end

[ntokens,tokens]=fltokens(fline);
itoken=1;
[tokens,itoken]=nextkey(fidi,tokens,itoken);
method=tokens{1}{itoken};
display(sprintf('Dakota method=%s.',method));

%%  loop through the file to find the Dakota variables

fseek(fidi,0,'bof');
[fline]=findline(fidi,'variables');
if ~ischar(fline)
    error('No Dakota variables in file %s.',filei);
end

[ntokens,tokens]=fltokens(fline);
itoken=1;
[dvar]=variables_parse(fidi,tokens,itoken);

%%  loop through the file to find the Dakota responses

fseek(fidi,0,'bof');
[fline]=findline(fidi,'responses');
if ~ischar(fline)
    error('No Dakota responses in file %s.',filei);
end

[ntokens,tokens]=fltokens(fline);
itoken=1;
[dresp]=responses_parse(fidi,tokens,itoken);

%%  loop through the file to find the Dakota response and probability levels
%   (even though they're in method section, process after responses)

fseek(fidi,0,'bof');
[fline]=findline(fidi,'method');

[ntokens,tokens]=fltokens(fline);
itoken=1;
[dresp]=resplevels(fidi,tokens,itoken,dresp);

%%  loop through the file to verify the end

display('End of file successfully reached.');
fclose(fidi);

end

%%  function to parse the dakota variables

function [dvar]=variables_parse(fidi,tokens,itoken)

display('Reading Dakota variables.');
dvar=[];
ncdv=0;
nnuv=0;
ncsv=0;

%  read next keyword

[tokens,itoken]=nextkey(fidi,tokens,itoken);
if ~itoken
    warning('variables_parse:empty',...
        'Dakota variables section is empty.');
end

%  process current keyword
%  (note that this is using dakota 4.1 keywords.  dakota 4.2
%  keywords are order-dependent.)

while itoken
    keyword=tokens{1}{itoken};
    display(sprintf('  Dakota keyword=%s.',keyword));

%  switch according to the keyword

    switch lower(keyword)
        case 'continuous_design'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            ncdv=tlist;
            dvar.cdv=[];
            display(sprintf('    Number of Dakota %s variables=%d.',...
                    'continuous_design',ncdv));
        case 'cdv_initial_point'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            for i=1:length(tlist)
                dvar.cdv(i).initpt    =tlist(i);
            end
        case 'cdv_lower_bounds'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            for i=1:length(tlist)
                dvar.cdv(i).lower     =tlist(i);
            end
        case 'cdv_upper_bounds'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            for i=1:length(tlist)
                dvar.cdv(i).upper     =tlist(i);
            end
        case 'cdv_descriptors'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            for i=1:length(tlist)
                dvar.cdv(i).descriptor=char(tlist(i));
            end

        case 'normal_uncertain'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            nnuv=tlist;
            dvar.nuv=[];
            display(sprintf('    Number of Dakota %s variables=%d.',...
                    'normal_uncertain',nnuv));
        case 'nuv_means'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            for i=1:length(tlist)
                dvar.nuv(i).mean      =tlist(i);
            end
        case 'nuv_std_deviations'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            for i=1:length(tlist)
                dvar.nuv(i).stddev    =tlist(i);
            end
        case 'nuv_lower_bounds'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            for i=1:length(tlist)
                dvar.nuv(i).lower     =tlist(i);
            end
        case 'nuv_upper_bounds'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            for i=1:length(tlist)
                dvar.nuv(i).upper     =tlist(i);
            end
        case 'nuv_descriptors'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            for i=1:length(tlist)
                dvar.nuv(i).descriptor=char(tlist(i));
            end

        case 'continuous_state'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            ncsv=tlist;
            dvar.csv=[];
            display(sprintf('    Number of Dakota %s variables=%d.',...
                    'continuous_state',ncsv));
        case 'csv_initial_state'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            for i=1:length(tlist)
                dvar.csv(i).initst    =tlist(i);
            end
        case 'csv_lower_bounds'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            for i=1:length(tlist)
                dvar.csv(i).lower     =tlist(i);
            end
        case 'csv_upper_bounds'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            for i=1:length(tlist)
                dvar.csv(i).upper     =tlist(i);
            end
        case 'csv_descriptors'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            for i=1:length(tlist)
                dvar.csv(i).descriptor=char(tlist(i));
            end

        otherwise
            warning('variables_parse:unrec_key',...
                'Unrecognized keyword ''%s''.',keyword);
            [tokens,itoken]=nextkey(fidi,tokens,itoken);
    end

%  check for eof or start of new section

    if (~itoken) || ...
       strncmpi(tokens{1}{itoken},'strategy' ,8) || ...
       strncmpi(tokens{1}{itoken},'method'   ,6) || ...
       strncmpi(tokens{1}{itoken},'model'    ,5) || ...
       strncmpi(tokens{1}{itoken},'variables',9) || ...
       strncmpi(tokens{1}{itoken},'interface',9) || ...
       strncmpi(tokens{1}{itoken},'responses',9)

%  supply default descriptors if necessary

        if isfield(dvar,'cdv') && ~isfield(dvar.cdv,'descriptor')
            for i=1:ncdv
                dvar.cdv(i).descriptor=sprintf('cdv_%d',i);
            end
        end
        if isfield(dvar,'nuv') && ~isfield(dvar.nuv,'descriptor')
            for i=1:nnuv
                dvar.nuv(i).descriptor=sprintf('nuv_%d',i);
            end
        end
        if isfield(dvar,'csv') && ~isfield(dvar.csv,'descriptor')
            for i=1:ncsv
                dvar.csv(i).descriptor=sprintf('csv_%d',i);
            end
        end
        return;
    end
end

end

%%  function to parse the dakota responses

function [dresp]=responses_parse(fidi,tokens,itoken)

display('Reading Dakota responses.');
dresp=[];
nof =0;
nlst=0;
nnic=0;
nnec=0;
nrf =0;

%  read next keyword

[tokens,itoken]=nextkey(fidi,tokens,itoken);
if ~itoken
    warning('responses_parse:empty',...
        'Dakota responses section is empty.');
end

%  process current keyword

while itoken
    keyword=tokens{1}{itoken};
    display(sprintf('  Dakota keyword=%s.',keyword));

%  switch according to the keyword

    switch lower(keyword)
        case 'num_objective_functions'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            nof =tlist;
            dresp.of =[];
            display(sprintf('    Number of Dakota %s=%d.',...
                    'objective_functions',nof));
        case 'objective_function_scale_types'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            for i=1:length(tlist)
                dresp.of(i).scale_type=char(tlist(i));
            end
        case 'objective_function_scales'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            for i=1:length(tlist)
                dresp.of(i).scale     =tlist(i);
            end
        case 'multi_objective_weights'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            for i=1:length(tlist)
                dresp.of(i).weight    =tlist(i);
            end

        case 'num_least_squares_terms'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            nlst=tlist;
            dresp.lst=[];
            display(sprintf('    Number of Dakota %s=%d.',...
                    'least_squares_terms',nlst));
        case 'least_squares_term_scale_types'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            for i=1:length(tlist)
                dresp.lst(i).scale_type=char(tlist(i));
            end
        case 'least_squares_term_scales'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            for i=1:length(tlist)
                dresp.lst(i).scale     =tlist(i);
            end
        case 'least_squares_weights'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            for i=1:length(tlist)
                dresp.lst(i).weight    =tlist(i);
            end

        case 'num_nonlinear_inequality_constraints'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            nnic=tlist;
            dresp.nic=[];
            display(sprintf('    Number of Dakota %s=%d.',...
                    'nonlinear_inequality_constraints',nnic));
        case 'nonlinear_inequality_scale_types'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            for i=1:length(tlist)
                dresp.nic(i).scale_type=char(tlist(i));
            end
        case 'nonlinear_inequality_scales'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            for i=1:length(tlist)
                dresp.nic(i).scale     =tlist(i);
            end
        case 'nonlinear_inequality_lower_bounds'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            for i=1:length(tlist)
                dresp.nic(i).lower     =tlist(i);
            end
        case 'nonlinear_inequality_upper_bounds'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            for i=1:length(tlist)
                dresp.nic(i).upper     =tlist(i);
            end

        case 'num_nonlinear_equality_constraints'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            nnec=tlist;
            dresp.nec=[];
            display(sprintf('    Number of Dakota %s=%d.',...
                    'nonlinear_equality_constraints',nnec));
        case 'nonlinear_equality_scale_types'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            for i=1:length(tlist)
                dresp.nec(i).scale_type=char(tlist(i));
            end
        case 'nonlinear_equality_scales'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            for i=1:length(tlist)
                dresp.nec(i).scale     =tlist(i);
            end
        case 'nonlinear_equality_targets'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            for i=1:length(tlist)
                dresp.nec(i).target    =tlist(i);
            end

        case 'num_response_functions'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            nrf =tlist;
            dresp.rf =[];
            display(sprintf('    Number of Dakota %s=%d.',...
                    'response_functions',nrf));

        case 'response_descriptors'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            desc=tlist;
        otherwise
            warning('responses_parse:unrec_key',...
                'Unrecognized keyword ''%s''.',keyword);
            [tokens,itoken]=nextkey(fidi,tokens,itoken);
    end

%  check for eof or start of new section

    if (~itoken) || ...
       strncmpi(tokens{1}{itoken},'strategy' ,8) || ...
       strncmpi(tokens{1}{itoken},'method'   ,6) || ...
       strncmpi(tokens{1}{itoken},'model'    ,5) || ...
       strncmpi(tokens{1}{itoken},'variables',9) || ...
       strncmpi(tokens{1}{itoken},'interface',9) || ...
       strncmpi(tokens{1}{itoken},'responses',9)

%  assign specified or supply default descriptors

        if exist('desc','var')
            idesc=0;
            if isfield(dresp,'of' )
                for i=1:nof
                    idesc=idesc+1;
                    dresp.of(i).descriptor=char(desc(idesc));
                end
            end
            if isfield(dresp,'lst')
                for i=1:nlst
                    idesc=idesc+1;
                    dresp.lst(i).descriptor=char(desc(idesc));
                end
            end
            if isfield(dresp,'nic')
                for i=1:nnic
                    idesc=idesc+1;
                    dresp.nic(i).descriptor=char(desc(idesc));
                end
            end
            if isfield(dresp,'nec')
                for i=1:nnec
                    idesc=idesc+1;
                    dresp.nec(i).descriptor=char(desc(idesc));
                end
            end
            if isfield(dresp,'rf' )
                for i=1:nrf
                    idesc=idesc+1;
                    dresp.rf(i).descriptor=char(desc(idesc));
                end
            end

        else
            if isfield(dresp,'of' )
                for i=1:nof
                    dresp.of(i).descriptor=sprintf('obj_fn_%d',i);
                end
            end
            if isfield(dresp,'lst')
                for i=1:nlst
                    dresp.lst(i).descriptor=sprintf('least_sq_term_%d',i);
                end
            end
            if isfield(dresp,'nic')
                for i=1:nnic
                    dresp.nic(i).descriptor=sprintf('nln_ineq_con_%d',i);
                end
            end
            if isfield(dresp,'nec')
                for i=1:nnec
                    dresp.nec(i).descriptor=sprintf('nln_eq_con_%d',i);
                end
            end
            if isfield(dresp,'rf' )
                for i=1:nrf
                    dresp.rf(i).descriptor=sprintf('response_fn_%d',i);
                end
            end
        end
        return;
    end
end

end

%%  function to read the number and levels of responses

function [dresp]=resplevels(fidi,tokens,itoken,dresp)

display('Reading Dakota response levels.');

%  read next keyword

[tokens,itoken]=nextkey(fidi,tokens,itoken);
if ~itoken
    warning('resplevels:empty',...
        'Dakota method section is empty.');
end

%  process current keyword

while itoken
    keyword=tokens{1}{itoken};
    display(sprintf('  Dakota keyword=%s.',keyword));

%  switch according to the keyword

    switch lower(keyword)
        case 'nond_sampling'
            [tokens,itoken]=nextkey(fidi,tokens,itoken);
        case 'nond_local_reliability'
            [tokens,itoken]=nextkey(fidi,tokens,itoken);
        case 'num_response_levels'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            nresp=tlist;
        case 'response_levels'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            nrespl=tlist;
        case 'num_probability_levels'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            nprob=tlist;
        case 'probability_levels'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            nprobl=tlist;
        case 'num_reliability_levels'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            nrel =tlist;
        case 'reliability_levels'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            nrell =tlist;
        case 'num_gen_reliability_levels'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            ngrel=tlist;
        case 'gen_reliability_levels'
            [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
            ngrell=tlist;
        case 'compute'
            [tokens,itoken]=nexttoken(fidi,tokens,itoken);
            [tokens,itoken]=nextkey(fidi,tokens,itoken);
        otherwise
            warning('resplevels:unrec_key',...
                'Unrecognized keyword ''%s''.',keyword);
            [tokens,itoken]=nextkey(fidi,tokens,itoken);
    end

%  check for eof or start of new section

    if (~itoken) || ...
       strncmpi(tokens{1}{itoken},'strategy' ,8) || ...
       strncmpi(tokens{1}{itoken},'method'   ,6) || ...
       strncmpi(tokens{1}{itoken},'model'    ,5) || ...
       strncmpi(tokens{1}{itoken},'variables',9) || ...
       strncmpi(tokens{1}{itoken},'interface',9) || ...
       strncmpi(tokens{1}{itoken},'responses',9)

%  assemble the lists by response

        if exist('nrespl','var') && isfield(dresp,'rf')
            if ~exist('nresp','var')
                nresp(1:length(dresp.rf))=floor(length(nrespl)/length(dresp.rf));
            end
            ilist=1;
            for i=1:length(dresp.rf)
                dresp.rf(i).respl=nrespl(ilist:ilist+nresp(i)-1);
                ilist=ilist+nresp(i);
            end
        end

        if exist('nprobl','var') && isfield(dresp,'rf')
            if ~exist('nprob','var')
                nprob(1:length(dresp.rf))=floor(length(nprobl)/length(dresp.rf));
            end
            ilist=1;
            for i=1:length(dresp.rf)
                dresp.rf(i).probl=nprobl(ilist:ilist+nprob(i)-1);
                ilist=ilist+nprob(i);
            end
        end

        if exist('nrell' ,'var') && isfield(dresp,'rf')
            if ~exist('nrel' ,'var')
                nrel (1:length(dresp.rf))=floor(length(nrell )/length(dresp.rf));
            end
            ilist=1;
            for i=1:length(dresp.rf)
                dresp.rf(i).rell =nrell (ilist:ilist+nrel (i)-1);
                ilist=ilist+nrel (i);
            end
        end

        if exist('ngrell','var') && isfield(dresp,'rf')
            if ~exist('ngrel','var')
                ngrel(1:length(dresp.rf))=floor(length(ngrell)/length(dresp.rf));
            end
            ilist=1;
            for i=1:length(dresp.rf)
                dresp.rf(i).grell=ngrell(ilist:ilist+ngrel(i)-1);
                ilist=ilist+ngrel(i);
            end
        end

        return;
    end
end

end

%%  function to find the next keyword

function [tokens,itoken]=nextkey(fidi,tokens,itoken)

%  start with next token

[tokens,itoken]=nexttoken(fidi,tokens,itoken);
if ~itoken
    return;
end

%  check for equal sign and skip subsequent list

if (itoken <= length(tokens{1})) && ...
   strncmp(tokens{1}{itoken},'=',1)
    [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken);
end

end

%%  function to find the next token

function [tokens,itoken]=nexttoken(fidi,tokens,itoken)

%  start with next token

itoken=itoken+1;

%  read next line if necessary

if (itoken > length(tokens{1}))
    fline=readline(fidi);
    if isempty(fline)
        tokens={};
        itoken=0;
        return;
    end
    [ntokens,tokens]=fltokens(fline);
    itoken=1;
end

end

%%  function to read a list of tokens

function [tlist,tokens,itoken]=readtlist(fidi,tokens,itoken)

%  start with next token (which should be equal sign, unless
%  equal sign was already read to determine existence of list)

itoken=itoken+1;

%  read next line if necessary

if (itoken > length(tokens{1}))
    fline=readline(fidi);
    if isempty(fline)
        tokens={};
        itoken=0;
        return;
    end
    [ntokens,tokens]=fltokens(fline);
    itoken=1;
end

%  check for equal sign and skip

if strncmp(tokens{1}{itoken},'=',1)
    itoken=itoken+1;
end

ilist=0;

%  accumulate list until non-numeric and non-quoted-string (or eof)
%  is encountered

while 1
    for i=itoken:length(tokens{1})
        if isnumeric(tokens{1}{i})
            ilist=ilist+1;
            tlist(ilist)=tokens{1}{i};
        elseif ischar(tokens{1}{i}) && ...
               (strncmp(tokens{1}{i}(1)  ,'''',1) && ...
                strncmp(tokens{1}{i}(end),'''',1)) || ...
               (strncmp(tokens{1}{i}(1)  ,'"',1) && ...
                strncmp(tokens{1}{i}(end),'"',1))
            ilist=ilist+1;
            tlist(ilist)=cellstr(tokens{1}{i}(2:end-1));
        else
            itoken=i;
            return
        end
    end
    fline=readline(fidi);
    if isempty(fline)
        tokens={};
        itoken=0;
        return;
    end
    [ntokens,tokens]=fltokens(fline);
    itoken=1;
end

end

%%  function to find a file line starting with a specified string

function [fline]=findline(fidi,string)

ipos=ftell(fidi);

while 1
    fline=readline(fidi);
    if isempty(fline)
        break;
    else
        if (strncmpi(fline,string,length(string)))
            return;
        end
    end
end

%  issue warning and reset file position

warning('findline:str_not_found',...
    'String ''%s'' not found in file.',string);
fseek(fidi,ipos,'bof');

end

%%  function to read a file line ignoring comments and blanks

function [fline]=readline(fidi)

while 1
    fline=fgetl(fidi);
    if ~ischar(fline)
        fline=[];
        return;
    end

    for ichar=1:length(fline)
        if ~strncmp(fline(ichar),' ',1) && ...
           ~strncmp(fline(ichar),'	',1)
            break;
        end
    end
    if isempty(fline) || ...
       (ichar > length(fline)) || ...
       strncmp(fline(ichar),'#',1)
        continue;
    else
        return;
    end
end

end

%%  function to parse a file line into tokens

function [ntokens,tokens]=fltokens(fline)

if ~ischar(fline)
    ntokens=-1;
    tokens={};
    return;
end
if isempty(fline)
    ntokens=0;
    tokens={};
    return;
end

strings=textscan(fline,'%s','delimiter',' :,');
%for i=1:length(strings{1})
%    display(sprintf('i=%d; strings{1}{%d}=%s',i,i,strings{1}{i}))
%end
ntokens=0;
tokens{1}{length(strings)}='';

for i=1:length(strings{1})
    if isempty(strings{1}{i})
        continue
    end
    ntokens=ntokens+1;
    inum=sscanf(strings{1}{i},'%f');
    if isempty(inum)
        tokens{1}{ntokens}=strings{1}{i};
%         display(sprintf('i=%d; tokens{1}{%d}=%s',...
%             i,ntokens,tokens{1}{ntokens}))
    else
        tokens{1}{ntokens}=inum;
%         display(sprintf('i=%d; tokens{1}{%d}=%f',...
%             i,ntokens,tokens{1}{ntokens}))
    end
end

end
