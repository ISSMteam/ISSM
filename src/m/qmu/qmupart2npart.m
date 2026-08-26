function npart=qmupart2npart(vector)
	%QMUPART2NPART - compute the number of partitions from a partition vector
	%
	%   vector is full of -1 (no partition) and 0 to npart. We need to identify npart.
	%
	%   Usage:
	%      npart=qmupart2npart(vector)
	npart=max(vector)+1;
