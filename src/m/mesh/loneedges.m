function edges=loneedges(md)

	edges=[];

	for e=1:md.mesh.numberofelements, 
		for j=1:3, 
			if j==3,
				ind1=md.mesh.elements(e,3);
				ind2=md.mesh.elements(e,1);
			else
				ind1=md.mesh.elements(e,j);
				ind2=md.mesh.elements(e,j+1);
			end

			%edge ind1 and ind2: 
			els1=md.mesh.vertexconnectivity(ind1,1: md.mesh.vertexconnectivity(ind1,end));
			els2=md.mesh.vertexconnectivity(ind2,1: md.mesh.vertexconnectivity(ind2,end));
			els=intersect(els1,els2);
			if length(els)~=2,
				edges=[edges;[ind1,ind2]];
			end
			end
		end
	end
