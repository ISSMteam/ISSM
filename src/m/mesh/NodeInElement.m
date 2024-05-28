function node_in_element=NodeInElement(newx,newy,elements,x,y,nodeconnectivity)
% NODEINELEMENT - find for a list of nodes (in newx,newy), which elements in the mesh (elements,x,y) they belong to.
%
%    Usage:
%      node_in_element=NodeInElement(newx,newy,elements,x,y,md.mesh.vertexconnectivity);
%
%  See also Nodeconnectivity
%
epsilon=10^-10;

%compute some quantities that will speed up the process
x3x1=x(elements(:,1))-x(elements(:,3));
y3y1=y(elements(:,1))-y(elements(:,3));
x3x2=x(elements(:,2))-x(elements(:,3));
y3y2=y(elements(:,2))-y(elements(:,3));
x3=x(elements(:,3));
y3=y(elements(:,3));
delta=x(elements(:,2)).*y(elements(:,3))-y(elements(:,2)).*x(elements(:,3))-x(elements(:,1)).*y(elements(:,3))+y(elements(:,1)).*x(elements(:,3))+x(elements(:,1)).*y(elements(:,2))-y(elements(:,1)).*x(elements(:,2));

%max connectivity:
max_connectivity=max(nodeconnectivity(:,end));
node_in_element=zeros(length(newx),max_connectivity+1); %last column is the number of elements to which the row node is connected.

for i=1:length(newx),
	x0=newx(i);
	y0=newy(i);

	%first area coordinate
	area_1=(y3y2.*(x0-x3)-x3x2.*(y0-y3))./delta;
	%second area coordinate
	area_2=(x3x1.*(y0-y3)-y3y1.*(x0-x3))./delta;
	%third area coordinate
	area_3=1-area_1-area_2;

	%get elements for which all area coordinates are positive (meaning (x0,y0) belongs to these elements
	pos=find((area_1>=0-epsilon) & (area_2>=0-epsilon) & (area_3>=0-epsilon));

	num_elements=length(pos);

	node_in_element(i,1:num_elements)=pos;
	node_in_element(i,end)=num_elements;
end
