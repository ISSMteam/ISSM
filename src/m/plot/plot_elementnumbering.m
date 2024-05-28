function plot_elementnumbering(md,options,width,i)
%PLOT_ELEMENTNUMBERING - plot element numbering
%
%   Usage:
%      plot_elementnumbering(md,options,width,i);
%
%   See also: PLOTMODEL, PLOT_UNIT, PLOT_MANAGER

subplot(width,width,i); 

%process data and model
[x y z elements is2d isplanet]=processmesh(md,[],options);
[elementnumbers datatype]=processdata(md,[1:md.mesh.numberofelements]',options);

%plot
if is2d
	%plot mesh 
	A=elements(:,1); B=elements(:,2); C=elements(:,3);
	patch( 'Faces', [A B C], 'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');

	%highlight
	pos=getfieldvalue(options,'highlight',[]);
	A=elements(pos,1); B=elements(pos,2); C=elements(pos,3);
	patch( 'Faces', [A B C], 'Vertices', [x y z],'FaceVertexCData', [0.9 0.5 0.5],'FaceColor','flat','EdgeColor','black');

	%numbering
	X   = sum(x(elements(:,1:3)),2)/3;
	Y   = sum(y(elements(:,1:3)),2)/3;
	Z   = sum(z(elements(:,1:3)),2)/3;
	NUM = elementnumbers;
	if exist(options,'axis')
		AXIS = getfieldvalue(options,'axis');
		pos = find(X>AXIS(1) & X<AXIS(2) & Y>AXIS(3) & Y<AXIS(4));
		X = X(pos); Y=Y(pos); Z=Z(pos); NUM=NUM(pos);
	end
	text(X,Y,Z,num2str(NUM),'HorizontalAlignment','center','VerticalAlignment','middle');
else
	%plot mesh 
	A=elements(:,1); B=elements(:,2); C=elements(:,3); D=elements(:,4); E=elements(:,5); F=elements(:,6);
	patch( 'Faces', [A B C],  'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
	patch( 'Faces', [D E F],  'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
	patch( 'Faces', [A B E D],'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
	patch( 'Faces', [B E F C],'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
	patch( 'Faces', [C A D F],'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');

	%highlight
	pos=getfieldvalue(options,'highlight',[]);
	A=elements(pos,1); B=elements(pos,2); C=elements(pos,3); D=elements(pos,4); E=elements(pos,5); F=elements(pos,6);
	patch( 'Faces', [A B C],  'Vertices', [x y z],'FaceVertexCData', [0.9 0.5 0.5],'FaceColor','flat','EdgeColor','black');
	patch( 'Faces', [D E F],  'Vertices', [x y z],'FaceVertexCData', [0.9 0.5 0.5],'FaceColor','flat','EdgeColor','black');
	patch( 'Faces', [A B E D],'Vertices', [x y z],'FaceVertexCData', [0.9 0.5 0.5],'FaceColor','flat','EdgeColor','black');
	patch( 'Faces', [B E F C],'Vertices', [x y z],'FaceVertexCData', [0.9 0.5 0.5],'FaceColor','flat','EdgeColor','black');
	patch( 'Faces', [C A D F],'Vertices', [x y z],'FaceVertexCData', [0.9 0.5 0.5],'FaceColor','flat','EdgeColor','black');

	%numbering
	X   = sum(x(elements(:,1:6)),2)/6;
	Y   = sum(y(elements(:,1:6)),2)/6;
	Z   = sum(z(elements(:,1:6)),2)/6;
	NUM = elementnumbers;
	if exist(options,'axis')
		AXIS = getfieldvalue(options,'axis');
		pos = find(X>AXIS(1) & X<AXIS(2) & Y>AXIS(3) & Y<AXIS(4));
		X = X(pos); Y=Y(pos); Z=Z(pos); NUM=NUM(pos);
	end
	text(X,Y,Z,num2str(NUM),'HorizontalAlignment','center','VerticalAlignment','middle');
end

%apply options
options=addfielddefault(options,'title','Element numbering');
options=addfielddefault(options,'colorbar',0);
applyoptions(md,[],options);
