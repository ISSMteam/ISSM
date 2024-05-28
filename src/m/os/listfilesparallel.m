function list=listfilesparallel(rank,numprocs)
%LISTFILESPARALLEL list files inside a directory, depending on rank  and number of processors running this routine.
%        this is very OS dependent.
%
%   usage: list=listfilesparallel(rank,numprocs);
%
%
%   see also LS DIR LISTFILES

list=listfiles';
numfiles=numel(list);

%we now have a list, split it between all the processors.
[i1,i2]=parallelrange(rank,numprocs,numfiles);
list=list(i1:i2);
