%
%  function to return the string dimension of an array element
%
%  function [sdim]=string_dim(a,idim,varargin)
%
function [sdim]=string_dim(a,idim,varargin)

if ~nargin
    help string_dim
    return
end

%  check for scalar

if (numel(a) == 1) && (idim == 1)
    sdim='';
    return
end

%  check for overflow

if (idim > numel(a))
    if ~isempty(inputname(1))
        error('Index %d exceeds number of elements in array ''%s''.',...
            idim,inputname(1));
    else
        error('Index %d exceeds number of elements in array %d.',...
            idim,1);
    end
end

%  check for column or row vector (Matlab uses a minimum of two
%  dimensions, so this won't match Matlab standard output)

for iarg=1:nargin-2
    if strcmpi(varargin{iarg},'vector')
        if (ndims(a) == 2) && ((size(a,1) == 1) || (size(a,2) == 1))
            sdim =['(' num2str(idim) ')'];
            return
        end
    end
end

%  do the general case

asize=size(a);
index=zeros(size(asize));
aprod=prod(asize);
idim =idim-1;

%  calculate indices base 0 and convert to base 1

%  note that ind2sub might be useful, except that it requires a list
%  of scalars rather than a vector for output.

for i=length(asize):-1:1
    aprod=aprod/asize(i);
    index(i)=floor(idim/aprod);
    idim=idim-index(i)*aprod;
end
index=index+1;

%  assemble string for output

sdim ='(';
for i=1:length(asize)-1;
    sdim =[sdim num2str(index(i)) ','];
end
sdim =[sdim num2str(index(end)) ')'];

end

