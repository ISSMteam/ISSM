function plot_edges(md,options,width,i,datai)
%PLOT_SEGMENTS - plot edges, with different colors according to segment markers.
%
%   Usage:
%      plot_edges(md,options,width,i);
%
%   See also: PLOTMODEL

%plot mesh boundaries
subplot(width,width,i); 

%process mesh and data
[x y z elements is2d isplanet]=processmesh(md,[],options);
edges=md.mesh.edges;
if isnan(edges)
	error('edges in NaN')
end

if dimension(md.mesh)==2,
	%plot mesh
	A=elements(:,1); B=elements(:,2); C=elements(:,3); 
	h1=patch( 'Faces', [A B C], 'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
	hold on;
	text(sum(x(edges(:,1:2)),2)/2,sum(y(edges(:,1:2)),2)/2,sum(z(edges(:,1:2)),2)/2,...
		num2str(transpose(1:size(edges,1))),...
		'backgroundcolor',[0.8 0.9 0.8],'HorizontalAlignment','center','VerticalAlignment','middle');
else
	error('plot_edges: 3d plot of edges not supported yet!');
end

%apply options
options=addfielddefault(options,'title','Edges');
options=addfielddefault(options,'colorbar',0);
applyoptions(md,[],options);
