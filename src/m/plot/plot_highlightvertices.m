function plot_highlightvertices(md,options,width,i)
%PLOT_HIGHLIGHTVERTICES - plot selected vertices
%
%   Usage:
%      plot_highlightvertices(md,options,width,i);
%
%   See also: PLOTMODEL

%process data and model
[x y z elements is2d isplanet]=processmesh(md,[],options);
[nodenumbers datatype]=processdata(md,[1:md.mesh.numberofvertices]',options);

%plot
subplot(width,width,i); 

if is2d
	%plot mesh 
	A=elements(:,1); B=elements(:,2); C=elements(:,3);
	patch( 'Faces', [A B C], 'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');

	%Highlight
	pos=getfieldvalue(options,'highlight',[]);
	text(x(pos),y(pos),z(pos),num2str(pos(:)),...
		'backgroundcolor',[1 0 0],'HorizontalAlignment','center','VerticalAlignment','middle');
else
	%plot mesh 
	A=elements(:,1); B=elements(:,2); C=elements(:,3); D=elements(:,4); E=elements(:,5); F=elements(:,6);
	patch( 'Faces', [A B C],  'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
	patch( 'Faces', [D E F],  'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
	patch( 'Faces', [A B E D],'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
	patch( 'Faces', [B E F C],'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
	patch( 'Faces', [C A D F],'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');

	%Highlight
	pos=getfieldvalue(options,'highlight',[]);
	text(x(pos),y(pos),z(pos),num2str(pos(:)),...
		'backgroundcolor',[1 0 0],'HorizontalAlignment','center','VerticalAlignment','middle');
end

%apply options
if ~exist(options,'highlight')
	disp('highlightvertices warning : highlight option empty, not node highlighted');
end
options=addfielddefault(options,'title','Highlighted Nodes');
options=addfielddefault(options,'colorbar',0);
applyoptions(md,[],options);
