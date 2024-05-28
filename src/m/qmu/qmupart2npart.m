function npart=qmupart2npart(vector)
	%vector is full of -1 (no partition) and 0 to npart. We need to identify npart=
	npart=max(vector)+1;
