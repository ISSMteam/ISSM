function factors=importancefactors(md,variablename,responsename,partition)
%IMPORTANCEFACTORS - compute importance factors for a certain variable and response.
%
%   Usage:
%      factors=importancefactors(md,variablename,responsename)
%
%
%   Example: factors=importancefactors(md,'drag','max_vel');
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
importancefactors=zeros(1,0);
count=0;
for i=1:nfun,
	desvar=responsefunctions.var{i};
	if strncmpi(desvar,variablename,variablenamelength),
		importancefactors(end+1)=responsefunctions.impfac(i);
		count=count+1;
	end
end

if count==0,
	error('importancefactors error message: either response does not exist, or importancefactors are empty');
end

if count==1, %we have scalar
	factors=importancefactors;
	return;
elseif count==max(partition+1)
	%distribute importance factor
	factors=importancefactors(partition'+1); %partition was created to index "c" style
else
	%distribute importance factor
	factors=importancefactors(partition'+1); %partition was created to index "c" style
end

%weight importancefactors by area
%if numel(factors)==md.mesh.numberofvertices,
%	%get areas for each vertex.
%	aire=GetAreas(md.mesh.elements,md.mesh.x,md.mesh.y);
%	num_elements_by_node=md.nodeconnectivity(:,end);
%	grid_aire=zeros(md.mesh.numberofvertices,1);
%	for i=1:md.mesh.numberofvertices,
%		for j=1:num_elements_by_node(i),
%			grid_aire(i)=grid_aire(i)+aire(md.nodeconnectivity(i,j));
%		end
%	end
%	factors=factors./grid_aire;
%end
