function new_flags=PropagateFlagsUntilDistance(md,flags,distance)
%PROPAGATEFLAGSUNTILDISTANCE
%
% Usage: 
%              flags=PropagateFlagsUntilDistance(md,flags,distance)
%
%

new_flags=flags;

%make 3d work in 2d: 
if dimension(md.mesh)==3,
	md.mesh.x=md.mesh.x2d;
	md.mesh.y=md.mesh.y2d;
	md.mesh.elements=md.mesh.elements2d;
end

%find elements that are at the border of flags: 
flag_elements=find(flags);
conn=md.mesh.elementconnectivity(flag_elements,:);
pos=find(conn);conn(pos)=~flags(conn(pos));
sum_conn=sum(conn,2);
border_elements=flag_elements(find(sum_conn>=1));

%average x and y over elements: 
x_elem=md.mesh.x(md.mesh.elements)*[1;1;1]/3;
y_elem=md.mesh.y(md.mesh.elements)*[1;1;1]/3;

while 1,

	%keep copy of new_flags for this loop: 
	new_flags_bak=new_flags;

	%extend new flags by connectivity
	pos=find(new_flags);

	connected_elements=md.mesh.elementconnectivity(pos,:);
	connected_elements=connected_elements(find(connected_elements));
	new_flags(connected_elements)=1;

	%get new elements: 
	new_elements=find(new_flags & ~new_flags_bak);
	if ~length(new_elements),
		%we are done!
		break;
	end

	%check which of these new elements are more than distance away from the border elements
	for i=1:length(new_elements),
		dist=sqrt(     (x_elem(border_elements)-x_elem(new_elements(i))).^2 + (y_elem(border_elements)-y_elem(new_elements(i))).^2)-distance;
		if ~any(dist<0)
			%none of the border elements are within distance, this element is outside out area of interest.
			%ensure this element never gets found again in the connectivity.
			pos=find(md.mesh.elementconnectivity==new_elements(i));
			md.mesh.elementconnectivity(pos)=0;
			%exclude this new element from the new_flags!
			new_flags(new_elements(i))=0;
		end
	end
end
