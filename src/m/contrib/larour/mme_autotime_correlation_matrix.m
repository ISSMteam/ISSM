function matrix=mme_autotime_correlation_matrix(mme,type)


	%Out of a multi model ensemble (nsamples x nsteps) of runs, build 
	%a temporal correlation matrix (of size nsteps x nsteps)

	nsamples=size(mme,1);
	nsteps=size(mme,2);

	%initialize with 1 in the diagonal: 
	matrix=eye(nsteps,nsteps);

	%go through time steps, and fill up the top part. 
	for i=1:nsteps,
		for j=i+1:nsteps,
			matrix(i,j)=corr(mme(:,i),mme(:,j),'Type',type);
			matrix(j,i)=matrix(i,j);
		end
	end

