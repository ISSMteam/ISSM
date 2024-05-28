function md=meshprocessoutsiderifts(md,domainoutline)
%MESHPROCESSOUTSIDERIFTS - process rifts when they touch the domain outline
%
%   Usage:
%      md=meshprocessoutsiderifts(md,domain)
%

%go through rifts, and figure out which ones touch the domain outline
for i=1:length(md.rifts.riftstruct),

	%first, flag nodes that belong to the domain outline
	flags=ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,domainoutline,'node',0);

	rift=md.rifts.riftstruct(i);
	tips=rift.tips;
	outsidetips=tips(find(flags(rift.tips)==0));

	%we have found outsidetips, tips that touch the domain outline. go through them
	for j=1:length(outsidetips),

		tip=outsidetips(j);
		%find tip in the segments, take first segment (there should be 2) that holds tip, 
		%and node_connected_to_tip is the other node on this segment:
		tipindex=find(rift.segments(:,1)==tip); 
		if length(tipindex),
			tipindex=tipindex(1);
			node_connected_to_tip=rift.segments(tipindex,2);
		else
			tipindex=find(rift.segments(:,2)==tip); tipindex=tipindex(1);
			node_connected_to_tip=rift.segments(tipindex,1);
		end

		%ok, we have the tip node, and the first node connected to it, on the rift. Now, 
		%identify all the elements that are connected to the tip, and that are on the same 
		%side of the rift.
		A=tip;
		B=node_connected_to_tip;

		elements=[];

		while flags(B), %as long as B does not belong to the domain outline, keep looking.
			%detect elements on edge A,B:
			edgeelements=ElementsFromEdge(md.mesh.elements,A,B);
			%rule out those we already detected
			already_detected=ismember(edgeelements,elements);
			nextelement=edgeelements(find(~already_detected));
			%add new detected element to the list of elements we are looking for.
			elements=[elements;nextelement];
			%new B:
			B=md.mesh.elements(nextelement,find(~ismember(md.mesh.elements(nextelement,:),[A B])));
		end

		%take the list of elements on one side of the rift that connect to the tip, 
		%and duplicate the tip on them, so as to open the rift to the outside.
		num=length(md.mesh.x)+1;
		md.mesh.x=[md.mesh.x;md.mesh.x(tip)];
		md.mesh.y=[md.mesh.y;md.mesh.y(tip)];
		md.mesh.numberofvertices=num;

		%replace tip in elements
		newelements=md.mesh.elements(elements,:);
		pos=find(newelements==tip);
		newelements(pos)=num;
		md.mesh.elements(elements,:)=newelements;
		md.rifts.riftstruct(i).tips=[md.rifts.riftstruct(i).tips num];

		%deal with segments
		tipsegments=find((md.mesh.segments(:,1)==tip) | (md.mesh.segments(:,2)==tip));
		for k=1:length(tipsegments),
			segment_index=tipsegments(k);
			pos=find(md.mesh.segments(segment_index,1:2)~=tip);
			other_node=md.mesh.segments(segment_index,pos);
			if ~isconnected(md.mesh.elements,other_node,tip),
				pos=find(md.mesh.segments(segment_index,1:2)==tip);
				md.mesh.segments(segment_index,pos)=num;
			end
		end
	end
end

%Fill in rest of fields:
md.mesh.numberofelements=length(md.mesh.elements);
md.mesh.numberofvertices=length(md.mesh.x);
md.mesh.vertexonboundary=zeros(length(md.mesh.x),1); md.mesh.vertexonboundary(md.mesh.segments(:,1:2))=1;
end

function flag=isconnected(elements,A,B)% {{{
	%ISCONNECTED: are two nodes connected by a triangulation?
	%
	%   Usage: flag=isconnected(elements,A,B)
	%
	%

	elements=ElementsFromEdge(elements,A,B);
	if isempty(elements),
		flag=0;
	else
		flag=1;
	end
end % }}}
