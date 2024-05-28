function md=setmask2(md,landname,floatingicename,groundedicename)
%GEOGRAPHY2 - establish land, ice sheet and ice shelf areas in a domains.
%
%   Usage:
%      md=setmask2(md,landname,floatingicename,groundedicename)
%
%   Examples:
%      md=setmask2(md,'LandName.exp','Iceshelves.exp','Islands.exp');

%Get assigned fields
x=md.mesh.x;
y=md.mesh.y;
elements=md.mesh.elements;

%recover elements and nodes on land.
if ischar(landname),
	[vertexonland,elementonland]=ContourToMesh(elements,x,y,landname,'element and node',2);
elseif isfloat(landname),
	if size(landname,1)~=md.mesh.numberofelements,
		error('Landname for area must be of same size as number of elements in model');
	end
	elementonland=landname;
	vertexonland=zeros(md.mesh.numberofvertices,1);
	vertexonland(md.mesh.elements(find(elementonland),:))=1;
else
	error('Invalid area option option');
end

%Now, build the connectivity tables for this mesh.
if size(md.mesh.vertexconnectivity,1)~=md.mesh.numberofvertices,
	md.mesh.vertexconnectivity=NodeConnectivity(md.mesh.elements,md.mesh.numberofvertices);
end
if size(md.mesh.elementconnectivity,1)~=md.mesh.numberofelements,
	md.mesh.elementconnectivity=ElementConnectivity(md.mesh.elements,md.mesh.vertexconnectivity);
end

%any element with 3 nodes on land should be on land:
elementsonwater=find(~elementonland);
wrongelements=elementsonwater(find(( vertexonland(md.mesh.elements(elementsonwater,1)) + vertexonland(md.mesh.elements(elementsonwater,2)) + vertexonland(md.mesh.elements(elementsonwater,3)) ...
                  )==3));
elementonland(wrongelements)=1;

%any element with its barycentre on land should be on land: (only if landname is an expfile)
if ischar(landname),
weights={[1;1;1],[2;1;1],[1;2;1],[1;1;2]};
	for i=1:length(weights),
		xelem=x(md.mesh.elements)*weights{i}/sum(weights{i});
		yelem=y(md.mesh.elements)*weights{i}/sum(weights{i});
	end
	baryonland=ContourToNodes(xelem,yelem,landname,1);
	pos=find(~baryonland); elementonland(pos)=0;
	pos=find(baryonland); elementonland(pos)=1;
end

%figure out which elements on land are actually in the middle of the ocean!
pos1=find(elementonland); 
connectedtoland=md.mesh.elementconnectivity(pos1,:);
pos=find(connectedtoland); connectedtoland(pos)=1-elementonland(connectedtoland(pos));
connectedtolandsum=sum(connectedtoland,2);
waterelements=pos1(find(connectedtolandsum==3));
elementonland(waterelements)=0;

%figure out which elements on water  are actually in the middle of the land!
pos1=find(~elementonland); 
connectedtowater=md.mesh.elementconnectivity(pos1,:);
pos=find(connectedtowater); connectedtowater(pos)=elementonland(connectedtowater(pos));
connectedtowatersum=sum(connectedtowater,2);
landelements=pos1(find(connectedtowatersum==3));
elementonland(landelements)=1;

%recover arrays of ice shelf nodes and elements, and ice sheet nodes and elements.
elementonfloatingice=FlagElements(md,floatingicename);
elementongroundedice=FlagElements(md,groundedicename);

%Because groundedice nodes and elements can be included into an floatingice, we need to update. Remember, all the previous 
%arrays come from domain outlines that can intersect one another: 
vertexonfloatingice=zeros(md.mesh.numberofvertices,1);
vertexongroundedice=zeros(md.mesh.numberofvertices,1);
elementonfloatingice=double((elementonfloatingice & ~elementongroundedice));
elementongroundedice=double(~elementonfloatingice);
vertexonfloatingice(md.mesh.elements(find(elementonfloatingice),:))=1;
vertexongroundedice(md.mesh.elements(find(elementongroundedice),:))=1;

%now correct, so that none of the floatingice and groundedice elements and nodes are in the water.
pos=find(~elementonland);
elementonfloatingice(pos)=0; 
elementongroundedice(pos)=0;

pos=find(~vertexonland);
vertexonfloatingice(pos)=0; 
vertexongroundedice(pos)=0;

%create vertexonwater and elementonwater: 
vertexonwater=double(~vertexonland);
elementonwater=double(~elementonland);

%correct for islands:
vertexonfloatingice=double(vertexonfloatingice & ~vertexongroundedice);
elementonfloatingice=double(elementonfloatingice & ~elementongroundedice);

%now, groundedices are everything except iceshelves and water
vertexongroundedice=double(~vertexonfloatingice & ~vertexonwater);
elementongroundedice=double(~elementonfloatingice & ~elementonwater);

%Deal with segments on neumann:

%Get current connectivity
mesh.elementconnectivity=md.mesh.elementconnectivity;

%put 0 for elements on water
pos=find(mesh.elementconnectivity);
mesh.elementconnectivity(pos)=mesh.elementconnectivity(pos).*(~elementonwater(mesh.elementconnectivity(pos)));

%put line of ones for elements on water
pos=find(elementonwater);
mesh.elementconnectivity(pos,:)=1;% line of ones for elements on water so they won't be considered

%resort lines (zeros must be at the last column for findsegments)
mesh.elementconnectivity=sort(mesh.elementconnectivity,2,'descend');

%call findsegments to build segment using THIS conectivity
md.mesh.segments=findsegments(md,'mesh.elementconnectivity',mesh.elementconnectivity);

%some final checks: 
%check that no node thinks it's on an ice shelf or ice sheet, and lies actually in the middle of the water.
nodesgrounded=find(~vertexonwater);
lengthconnectivity=size(md.mesh.vertexconnectivity,2);
groundedcounters=md.mesh.vertexconnectivity(nodesgrounded,lengthconnectivity);
groundedconnectivity=md.mesh.vertexconnectivity(nodesgrounded,1:lengthconnectivity-1);
pos=find(groundedconnectivity);
groundedconnectivity(pos)=elementonwater(groundedconnectivity(pos));
groundedsum=sum(groundedconnectivity,2);
errorflags=find(groundedsum==groundedcounters);
errornodes=nodesgrounded(errorflags);

vertexonwater(errornodes)=1;
vertexongroundedice(errornodes)=0;
vertexonfloatingice(errornodes)=0;

%Return: 
md.mesh.segmentmarkers(:)=1;
