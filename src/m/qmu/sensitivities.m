function sens=sensitivies(md,variablename,responsename)
%SENSITIVIES - compute sensitivities for a certain variable and response.
%
%   Usage:
%      sens=sensitivities(md,variablename,responsename)
%
%
%   Example: sens=sensitivities(md,'DragCoefficient','MaxVel');
%

variablenamelength=length(variablename);

%go through all response functions and find the one corresponding to the correct responsename
responsefunctions=md.qmu.results.dresp_out;
found=0;
for i=1:length(responsefunctions),
	if strcmpi(responsefunctions(i).descriptor,responsename),
		found=i;
		break;
	end
end
if ~found,
	error('importancefactors error message: could not find correct response function');
end
responsefunctions=responsefunctions(found);
nfun=size(responsefunctions.var,1);

%Now recover response to the correct design variable
rawsens=zeros(0,1);
count=0;
for i=1:nfun,
	desvar=responsefunctions.var{i};
	if strncmpi(desvar,variablename,variablenamelength),
		rawsens(end+1,1)=responsefunctions.sens(i);
		count=count+1;
	end
end

%Now, if this was a distributed variable, the sensitivities need to be scaled by means of the input variable.
if IsScaled(variablename),

	%ipick up the variable in the model
	switch variablename,
		case 'thickness', variable = md.geometry.thickness; 
		otherwise, error(['scaled variable ' variablename  ' not associated to any model field']);
	end

	%average it onto the partition
	average_variable=AreaAverageOntoPartition(md,variable);

	%scale the sensitivities: only where the average_variable is not 0 
	if ~isempty(rawsens),
		pos=find(average_variable);
		rawsens(pos)=rawsens(pos)./average_variable(pos);
	end
end

if count==0,
	error('sensitivities error message: either response does not exist, or sensitivities are empty');
end

if count==1, %we have scalar
	sens=rawsens;
	return;
else
	%project the sensitivities from the partition onto the mesh
	sens=rawsens(md.qmu.partition'+1); %md.qmu.partition was created to index "c" style
end
