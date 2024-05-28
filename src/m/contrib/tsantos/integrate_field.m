function total=integrate_field(index,x,y,field)
%INTEGRATE_FIELD: integrate a field over a 2D mesh
%
%   Usage:
%      total=integrate_field(index,x,y,field);
%
%   Examples:
%      volume=integrate_field(md.mesh.elements,md.mesh.x,md.mesh.y,md.geometry.thickness);

% areas of each element
x1=x(index(:,1)); x2=x(index(:,2)); x3=x(index(:,3));
y1=y(index(:,1)); y2=y(index(:,2)); y3=y(index(:,3));
areas=(0.5*((x2-x1).*(y3-y1)-(y2-y1).*(x3-x1)));

% element-wise integration
v1=index(:,1);	v2=index(:,2);	v3=index(:,3);
elem_int=areas.*mean(field([v1 v2 v3]),2);

% compute integration
total=sum(elem_int);

end
