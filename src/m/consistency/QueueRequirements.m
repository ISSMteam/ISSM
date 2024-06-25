function QueueRequirements(available_queues,queue_requirements_time,queue_requirements_np,queue,np,time)
%QUEUEREQUIREMENTS - queue requirements in time, number of cpus, by name of queue.
%
%   Usage: 
%      QueueRequirements(available_queues,queue_requirements_time,queue_requirements_np,np,time)

%Ok, go through requirements for current queue:
index=ismemberi(queue,available_queues);
if  ~index,
	%ok, either we a generic cluster, with 'none' queue, or we could not find the queue requirements
	if strcmpi(available_queues{1},'none'),
		%reset index to 1, so we can fish the requirements
		index=1;
	else
		string=available_queues{1};
		for i=2:length(available_queues),
			string=[string ' ' available_queues{i}];
		end
		error(['QueueRequirements error message: available queues are ' string]);
	end
end

%check on time requirements
rtime=queue_requirements_time(index);
if time<=0,
	error('QueueRequirements: time should be a positive number');
end
if time>rtime,
	error(['QueueRequirements: time should be < ' num2str(rtime) ' for queue: ' queue]);
end

%check on np requirements
if np<=0,
	error('QueueRequirements: np should be a positive number');
end
