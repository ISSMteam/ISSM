function mesh_profile=ProfileProjectOntoMesh(md,profile)
%PROFILEPROJECTONTOMESH: project a profile (made of arbitrary points) onto a mesh, so that we end 
%                        up with a list of segments self contained onto elements.
%
% Usage: mesh_profile=ProfileProjectOntoMesh(md,profile)
%
% See also intersections.m

%make a curve out of the mesh, to use the intersections routine.
rows=[md.mesh.elements md.mesh.elements(:,1)]'; rows=rows(:);
x=md.mesh.x(rows);
y=md.mesh.y(rows);

%[x0,y0] = intersections(profile.x,profile.y,x,y,1);
[x0,y0,indices,j] = intersections(profile.x,profile.y,x,y);

%  sort intersections to create segments in order and continuous along profile
[indices,isort]=sort(indices);
j =j (isort);
x0=x0(isort);
y0=y0(isort);

%process x0,y0 so they do not include profile.x or profile.y
processed_indices=[];
processed_x=[];
processed_y=[];
for i=1:numel(indices),
	if(((indices(i)-floor(indices(i)))~=0) && ((ceil(indices(i))-indices(i))~=0))
		processed_indices=[processed_indices;floor(indices(i))];
		processed_x=[processed_x;x0(i)];
		processed_y=[processed_y;y0(i)];
	end
end

%now merge profile.x,profile.y with processed_x,processed_y, at locations processed_indices:
newx=profile.x;
newy=profile.y;

count=1;
for i=1:numel(profile.x),
	pos=find(processed_indices==i);
	if ~isempty(pos),
		newx=[newx(1:count); processed_x(pos); newx(count+1:end)];
		newy=[newy(1:count); processed_y(pos); newy(count+1:end)];
		count=count+length(pos)+1;
	end
end

%now, for each node, figure out which element it belongs to.
node_in_element=NodeInElement(newx,newy,md.mesh.elements,md.mesh.x,md.mesh.y,md.mesh.vertexconnectivity);

% eliminate nodes that don't fall in any element
% (profile may start and/or end externally and/or cross holes in the model)

ind=find(node_in_element(:,end)>0);
newx=newx(ind,:);
newy=newy(ind,:);
node_in_element=node_in_element(ind,:);

mesh_profile=[newx(1:end-1) newy(1:end-1) newx(2:end) newy(2:end) zeros(length(newy(2:end)),1)];

%find which element each segment belongs to.
for i=1:length(newx)-1,
	common=intersect(node_in_element(i,1:node_in_element(i,end)), node_in_element(i+1,1:node_in_element(i+1,end)));
	mesh_profile(i,end)=common(1);
end
