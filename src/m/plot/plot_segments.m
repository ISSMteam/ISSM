function plot_segments(md,options,width,i,datai)
%PLOT_SEGMENTS - plot segments, with different colors according to segment markers.
%
%   Usage:
%      plot_segments(md,options,width,i);
%
%   See also: PLOTMODEL

%plot mesh boundaries
subplot(width,width,i); 

%process mesh and data
[x y z elements is2d isplanet]=processmesh(md,[],options);
segments=md.mesh.segments;

if dimension(md.mesh)==2,
	%plot mesh
	A=elements(:,1); B=elements(:,2); C=elements(:,3); 
	h1=patch( 'Faces', [A B C], 'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
	hold on;

	%highlight elements on neumann
	pos=segments(:,end);
	A=elements(pos,1); B=elements(pos,2); C=elements(pos,3); 
	h2=patch( 'Faces', [A B C], 'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','green','EdgeColor','black');

	if strcmpi(getfieldvalue(options,'segmentnumbering','off'),'on'),
		text(sum(x(segments(:,1:2)),2)/2,sum(y(segments(:,1:2)),2)/2,sum(z(segments(:,1:2)),2)/2+1,...
			num2str(md.mesh.segmentmarkers),...
			'backgroundcolor',[0.8 0.9 0.8],'HorizontalAlignment','center','VerticalAlignment','middle');
	end

	%display arrows pointing outward
	xstart=mean(x(segments(:,1:end-1)),2);
	ystart=mean(y(segments(:,1:end-1)),2);
	length=sqrt((x(segments(:,1))-x(segments(:,2))).^2 + (y(segments(:,1))-y(segments(:,2))).^2 );
	normal(:,1)=cos(atan2((x(segments(:,1))-x(segments(:,2))) , (y(segments(:,2))-y(segments(:,1)))));
	normal(:,2)=sin(atan2((x(segments(:,1))-x(segments(:,2))) , (y(segments(:,2))-y(segments(:,1)))));
	xend=xstart+length.*normal(:,1);
	yend=ystart+length.*normal(:,2);
	q=quiver(xstart,ystart,xend-xstart,yend-ystart); hold on;
	h3=plot(xstart,ystart,'r*');

else
	error('plot_segments: 3d plot of segments not supported yet!');
end

%legend (disable warnings)
warning off
legend([h2,q],'element on segment','normal vectors')
warning on

%apply options
options=addfielddefault(options,'title','Segment boundaries');
options=addfielddefault(options,'colorbar',0);
applyoptions(md,[],options);
