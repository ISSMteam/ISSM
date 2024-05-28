%
%  function to write a vector on multiple lines
%
function []=vector_write(fidi,sbeg,vec,nmax,cmax)

if ~exist('nmax','var') || isempty(nmax)
    nmax=Inf;
end
if ~exist('cmax','var') || isempty(cmax)
    cmax=Inf;
end

%  set up first iteration

svec =[];
nitem=nmax;
lsvec=cmax;

%  transpose vector from column-wise to row-wise

vec=vec';

%  assemble each line, flushing when necessary
for i=1:numel(vec)
    if isnumeric(vec(i))
        sitem=sprintf('%g'    ,vec(i));
    else
        sitem=sprintf('''%s''',char(vec(i)));
    end
    nitem=nitem+1;
    lsvec=lsvec+1+length(sitem);

    if (nitem <= nmax) && (lsvec <= cmax)
        svec=[svec ' ' sitem];
    else
        if ~isempty(svec)
            fprintf(fidi,'%s\n',svec);
        end
        svec=[sbeg sitem];
        nitem=1;
        lsvec=length(svec);
    end
end

%  flush buffer at end, if necessary

if ~isempty(svec)
    fprintf(fidi,'%s\n',svec);
end

end
