function plot_mesh(md,options,nlines,ncols,i)
%PLOT_MESH - plot model mesh
%
%   Usage:
%      plot_mesh(md,options,nlines,ncols,i);
%
%   See also: PLOTMODEL

%process data and model
[x y z elements is2d isplanet]=processmesh(md,[],options);

%plot mesh
subplot(nlines,ncols,i); 

%retrieve some options
linewidth=getfieldvalue(options,'linewidth',1);
edgecolor=getfieldvalue(options,'edgecolor','black');

%plot mesh
if is2d
	A=elements(:,1); B=elements(:,2); C=elements(:,3); 
	patch( 'Faces', [A B C], 'Vertices', [x y z],'FaceVertexCData',zeros(size(x)),'FaceColor','none','EdgeColor',edgecolor,'linewidth',linewidth);
else
	if size(elements,2)==6, %prisms
		A=elements(:,1); B=elements(:,2); C=elements(:,3); D=elements(:,4); E=elements(:,5); F=elements(:,6);
		patch( 'Faces', [A B C],  'Vertices', [x y z],'FaceVertexCData',zeros(size(x)),'FaceColor','none','EdgeColor',edgecolor,'linewidth',linewidth);
		patch( 'Faces', [D E F],  'Vertices', [x y z],'FaceVertexCData',zeros(size(x)),'FaceColor','none','EdgeColor',edgecolor,'linewidth',linewidth);
		patch( 'Faces', [A B E D],'Vertices', [x y z],'FaceVertexCData',zeros(size(x)),'FaceColor','none','EdgeColor',edgecolor,'linewidth',linewidth);
		patch( 'Faces', [B E F C],'Vertices', [x y z],'FaceVertexCData',zeros(size(x)),'FaceColor','none','EdgeColor',edgecolor,'linewidth',linewidth);
		patch( 'Faces', [C A D F],'Vertices', [x y z],'FaceVertexCData',zeros(size(x)),'FaceColor','none','EdgeColor',edgecolor,'linewidth',linewidth);
	elseif size(elements,2)==4, %tetras
		A=elements(:,1); B=elements(:,2); C=elements(:,3); D=elements(:,4);
		patch( 'Faces',[A B C],'Vertices', [x y z],'FaceVertexCData',zeros(size(x)),'FaceColor','none','EdgeColor',edgecolor,'linewidth',linewidth);
		patch( 'Faces',[A B D],'Vertices', [x y z],'FaceVertexCData',zeros(size(x)),'FaceColor','none','EdgeColor',edgecolor,'linewidth',linewidth);
		patch( 'Faces',[B C D],'Vertices', [x y z],'FaceVertexCData',zeros(size(x)),'FaceColor','none','EdgeColor',edgecolor,'linewidth',linewidth);
		patch( 'Faces',[C A D],'Vertices', [x y z],'FaceVertexCData',zeros(size(x)),'FaceColor','none','EdgeColor',edgecolor,'linewidth',linewidth);
	else %triangles (planet)
		A=elements(:,1); B=elements(:,2); C=elements(:,3); 
		if (size(elements,2)==4), D=elements(:,4); else D=C; end
		patch( 'Faces', [A B C D],  'Vertices', [x y z],'FaceVertexCData',zeros(size(x)),'FaceColor','none','EdgeColor',edgecolor,'linewidth',linewidth);
	end
end

%apply options
options=addfielddefault(options,'title','Mesh');
options=addfielddefault(options,'colorbar',0);
applyoptions(md,[],options);
