function plot_rifts(md,options,nlines,ncols,index)
%PLOT_RIFTS - plot rifts in a mesh
%
%   Usage:
%      plot_rifts(md,options,width,i);
%
%   See also: PLOTMODEL

%process data and model
[x y z elements is2d isplanet]=processmesh(md,[],options);

%plot mesh
subplot(nlines,ncols,index); 

%offset to separate rift flanks.
offset=getfieldvalue(options,'offset',500);
if isstruct(md.rifts.riftstruct),

	for i=1:size(md.rifts.riftstruct,1),
		penaltypairs=md.rifts.riftstruct(i).penaltypairs;

		normal=zeros(2,1);
		for j=1:size(penaltypairs,1),
			normal(1)=penaltypairs(j,5);
			normal(2)=penaltypairs(j,6);
			x(penaltypairs(j,1))=x(penaltypairs(j,1))-normal(1)*offset;
			y(penaltypairs(j,1))=y(penaltypairs(j,1))-normal(2)*offset;
		end
		if length(md.rifts.riftstruct(i).tips)==3,
			tip=md.rifts.riftstruct(i).tips(3);
			%who is tip connected to? 
			if isconnected(md.mesh.elements,penaltypairs(1,1),tip),
				normal(1)=penaltypairs(1,5);
				normal(2)=penaltypairs(1,6);
				x(tip)=x(tip)-normal(1)*offset;
				y(tip)=y(tip)-normal(2)*offset;
			end

			if isconnected(md.mesh.elements,penaltypairs(1,2),tip),
				normal(1)=penaltypairs(1,5);
				normal(2)=penaltypairs(1,6);
				x(tip)=x(tip)+normal(1)*offset;
				y(tip)=y(tip)+normal(2)*offset;
			end
			if isconnected(md.mesh.elements,penaltypairs(end,1),tip),
				normal(1)=penaltypairs(end,5);
				normal(2)=penaltypairs(end,6);
				x(tip)=x(tip)-normal(1)*offset;
				y(tip)=y(tip)-normal(2)*offset;
			end
			if isconnected(md.mesh.elements,penaltypairs(end,2),tip),
				normal(1)=penaltypairs(end,5);
				normal(2)=penaltypairs(end,6);
				x(tip)=x(tip)+normal(1)*offset;
				y(tip)=y(tip)+normal(2)*offset;
			end
		end
	end
end

%plot mesh
if is2d
	A=elements(:,1); B=elements(:,2); C=elements(:,3); 
	patch( 'Faces', [A B C], 'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
elseif isplanet,
	A=elements(:,1); B=elements(:,2); C=elements(:,3); 
	patch( 'Faces', [A B C], 'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
else
	A=elements(:,1); B=elements(:,2); C=elements(:,3); D=elements(:,4); E=elements(:,5); F=elements(:,6);
	patch( 'Faces', [A B C], 'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
	patch( 'Faces', [D E F], 'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
	patch( 'Faces', [A B E D], 'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
	patch( 'Faces', [B E F C ], 'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
	patch( 'Faces', [C A D F ], 'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
end

%apply options
options=addfielddefault(options,'title','Rifts');
options=addfielddefault(options,'colorbar',0);
applyoptions(md,[],options);
