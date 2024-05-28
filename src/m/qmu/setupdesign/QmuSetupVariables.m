function dvar=QmuSetupVariables(md,dvar,variables)

%get descriptor
descriptor=variables.descriptor;

%decide whether this is a distributed variable, which will drive whether we expand it into npart values,
%or if we just carry it forward as is. 

%ok, key off according to type of descriptor:
if strncmp(descriptor,'scaled_',7),
	%we have a scaled variable, expand it over the partition. First recover the partition: 
	partition=variables.partition;
	%figure out number of partitions: 
	npart=qmupart2npart(partition);
	%figure out number of time steps: 
	nt=variables.nsteps;

	if isa(variables,'uniform_uncertain'),
		nlower=size(variables.lower,1);
		nupper=size(variables.upper,1);
		if (nlower ~= npart || nupper ~=npart),
			error('QmuSetupVariables error message: upper and lower fields should have the same number of rows as the number of partitions');
		end
		nlower=size(variables.lower,2);
		nupper=size(variables.upper,2);
		if (nlower ~= nt || nupper ~= nt),
			error('QmuSetupVariables error message: upper and lower fields should have the same number of cols as the number of time steps');
		end
	elseif isa(variables,'normal_uncertain'),
		nstddev=size(variables.stddev,1);
		nmean=size(variables.mean,1);
		if (nstddev ~= npart || nmean ~=npart),
			error('QmuSetupVariables error message: stddev and mean fields should have the same number of rows as the number of partitions');
		end
		nstddev=size(variables.stddev,2);
		nmean=size(variables.mean,2);
		if (nstddev ~= nt || nmean ~=nt),
			error('QmuSetupVariables error message: stddev and mean fields should have the same number of cols as the number of time steps');
		end

	end

	%ok, dealing with semi-discrete distributed variable. Distribute according to how many 
	%partitions we want, and number of time steps: 
	if nt==1,
		for j=1:npart,
			dvar(end+1)           =variables;
			dvar(end  ).descriptor=sprintf('%s_%d',variables.descriptor,j);
			if isa(variables,'uniform_uncertain'),
				dvar(end  ).lower=variables.lower(j);
				dvar(end  ).upper=variables.upper(j);
			elseif isa(variables,'normal_uncertain'),
				dvar(end  ).stddev=variables.stddev(j);
				dvar(end  ).mean=variables.mean(j);
			end
		end
	else
		for j=1:npart,
			for k=1:nt,
				dvar(end+1)           =variables;
				dvar(end  ).descriptor=sprintf('%s_%d_%d',variables.descriptor,j,k);
				if isa(variables,'uniform_uncertain'),
					dvar(end  ).lower=variables.lower(j,k);
					dvar(end  ).upper=variables.upper(j,k);
				elseif isa(variables,'normal_uncertain'),
					dvar(end  ).stddev=variables.stddev(j,k);
					dvar(end  ).mean=variables.mean(j,k);
				end
			end
		end

	end

elseif strncmp(descriptor,'distributed_',12),

	%we have a distributed variable, expand it over the partition. First recover the partition: 
	partition=variables.partition;
	%figure out number of partitions: 
	npart=qmupart2npart(partition);

	if isa(variables,'uniform_uncertain'),
		nlower=size(variables.lower,1);
		nupper=size(variables.upper,1);
		if (nlower ~= npart || nupper ~=npart),
			error('QmuSetupVariables error message: upper and lower fields should have the same number of rows as the number of partitions');
		end
	elseif isa(variables,'normal_uncertain'),
		nstddev=size(variables.stddev,1);
		nmean=size(variables.mean,1);
		if (nstddev ~= npart || nmean ~=npart),
			error('QmuSetupVariables error message: stddev and mean fields should have the same number of rows as the number of partitions');
		end
	elseif isa(variables,'histogram_bin_uncertain'),
		ncounts=length(variables.counts);
		npairs=length(variables.pairs_per_variable);
		nabs=length(variables.abscissas);
		if (ncounts ~= npart),
			error(sprintf('QmuSetupVariables error message: counts size (%i) should be equal to the number of partitions (%i)',ncounts,npart));
		end
		if (npairs ~= npart),
			error(sprintf('QmuSetupVariables error message: pairs_per_variable size (%i) should be equal to the number of partitions (%i)',npairs,npart));
		end
		if (nabs ~= npart),
			error(sprintf('QmuSetupVariables error message: abscissas size (%i) should be equal to the number of partitions (%i)',nabs,npart));
		end
	end

	%ok, dealing with distributed variable. Distribute according to how many 
	%partitions we want.
	for j=1:npart,
		dvar(end+1)           =variables;
		dvar(end).descriptor=sprintf('%s_%d',variables.descriptor,j);
		if isa(variables,'uniform_uncertain'),
			dvar(end  ).lower=variables.lower(j);
			dvar(end  ).upper=variables.upper(j);
			dvar(end  ).partition=[];
		elseif isa(variables,'normal_uncertain'),
			dvar(end  ).stddev=variables.stddev(j);
			dvar(end  ).mean=variables.mean(j);
			dvar(end  ).partition=[];
		elseif isa(variables,'histogram_bin_uncertain'),
			dvar(end).pairs_per_variable=variables.pairs_per_variable(j);
			dvar(end).abscissas=variables.abscissas{j};
			dvar(end).counts=variables.counts{j};
			dvar(end).partition=[];
		end
	end
else
	dvar(end+1)=variables;
end
