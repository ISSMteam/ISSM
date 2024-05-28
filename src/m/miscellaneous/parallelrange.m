function [i1,i2]=parallelrange(rank,numprocs,globalsize)
%PARALLELRANGE - from a rank, and a number of processors, figure out a range, for parallel tasks.
%
%   Usage: 
%      [i1,i1]=parallelrange(rank,numprocs,globalsize)

num_local_rows=zeros(numprocs,1);

for i=1:numprocs,
	%we use floor. we under distribute rows. The rows left  are then redistributed, therefore resulting in a more even distribution.
	num_local_rows(i)=floor(globalsize/numprocs);
end

%There may be some rows left. Distribute evenly.
row_rest=globalsize - numprocs*floor(globalsize/numprocs);

for i=1:row_rest,
	num_local_rows(i)=num_local_rows(i)+1;
end

i1=sum(num_local_rows(1:rank-1))+1;
i2=i1+num_local_rows(rank)-1;
