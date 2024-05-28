function dresp=QmuSetupResponses(md,dresp,responses)

%get descriptor
descriptor=responses.descriptor;

%decide whether this is a distributed response, which will drive whether we expand it into npart values,
%or if we just carry it forward as is. 

%ok, key off according to type of descriptor:
if strncmp(descriptor,'scaled_',7),
	%we have a scaled response, expand it over the partition.

	%ok, dealing with semi-discrete distributed response. Distribute according to how many 
	%partitions we want
	npart=qmupart2npart(responses.partition);

	for j=1:npart,
		dresp(end+1)           =responses;
		dresp(end  ).descriptor=sprintf('%s_%d',responses.descriptor,j);
	end

else
	dresp(end+1)=responses;
end
