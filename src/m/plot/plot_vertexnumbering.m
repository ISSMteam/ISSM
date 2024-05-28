function plot_vertexnumbering(md,options,width,i)
%PLOT_VERTEXNUMBERING - plot vertex numbering
%
%   Usage:
%      plot_vertexnumbering(md,options,width,i);
%
%   See also: PLOTMODEL

%process data and model
[x y z elements is2d isplanet]=processmesh(md,[],options);
[vertexnumbers datatype]=processdata(md,[1:md.mesh.numberofvertices]',options);

%plot
subplot(width,width,i); 

if is2d
	%plot mesh 
	A=elements(:,1); B=elements(:,2); C=elements(:,3);
	patch( 'Faces', [A B C], 'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');

	%numbering
	X   = x;
	Y   = y;
	NUM = vertexnumbers;
	if exist(options,'axis')
		AXIS = getfieldvalue(options,'axis');
		pos = find(X>AXIS(1) & X<AXIS(2) & Y>AXIS(3) & Y<AXIS(4));
		X = X(pos); Y=Y(pos); NUM=NUM(pos);
	end
	text(X,Y,num2str(NUM),'backgroundcolor',[0.8 0.9 0.8],'HorizontalAlignment','center','VerticalAlignment','middle');

	%Highlight
	pos=getfieldvalue(options,'highlight',[]);
	text(x(pos),y(pos),z(pos),num2str(pos(:)),...
		'backgroundcolor',[1 0 0],'HorizontalAlignment','center','VerticalAlignment','middle');
else
	if size(elements,2)==6, %prisms
		%plot mesh 
		A=elements(:,1); B=elements(:,2); C=elements(:,3); D=elements(:,4); E=elements(:,5); F=elements(:,6);
		patch( 'Faces', [A B C],  'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
		patch( 'Faces', [D E F],  'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
		patch( 'Faces', [A B E D],'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
		patch( 'Faces', [B E F C],'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
		patch( 'Faces', [C A D F],'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
	elseif size(elements,2)==4, %tetras
		A=elements(:,1); B=elements(:,2); C=elements(:,3); D=elements(:,4);
		patch( 'Faces',[A B C],'Vertices', [x y z],'FaceVertexCData',zeros(size(x)),'FaceColor','none','EdgeColor','black');
		patch( 'Faces',[A B D],'Vertices', [x y z],'FaceVertexCData',zeros(size(x)),'FaceColor','none','EdgeColor','black');
		patch( 'Faces',[B C D],'Vertices', [x y z],'FaceVertexCData',zeros(size(x)),'FaceColor','none','EdgeColor','black');
		patch( 'Faces',[C A D],'Vertices', [x y z],'FaceVertexCData',zeros(size(x)),'FaceColor','none','EdgeColor','black');
	end

	%numbering
	X   = x;
	Y   = y;
	Z   = z;
	NUM = vertexnumbers;
	if exist(options,'axis')
		AXIS = getfieldvalue(options,'axis');
		pos = find(X>AXIS(1) & X<AXIS(2) & Y>AXIS(3) & Y<AXIS(4));
		X = X(pos); Y=Y(pos); Z=Z(pos); NUM=NUM(pos);
	end
	text(X,Y,Z,num2str(NUM),'backgroundcolor',[0.8 0.9 0.8],'HorizontalAlignment','center','VerticalAlignment','middle');

	%Highlight
	pos=getfieldvalue(options,'highlight',[]);
	text(x(pos),y(pos),z(pos),num2str(pos(:)),...
		'backgroundcolor',[1 0 0],'HorizontalAlignment','center','VerticalAlignment','middle');
end

%apply options
options=addfielddefault(options,'title','Node numbering');
options=addfielddefault(options,'colorbar',0);
applyoptions(md,[],options);
