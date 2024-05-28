function plot_highlightelements(md,options,width,i)
%PLOT_HIGHLIGHTELEMENTS - plot selected elements
%
%   Usage:
%      plot_highlightelements(md,options,width,i);
%
%   See also: PLOTMODEL

%plot mesh boundaries
subplot(width,width,i); 

%process data and model
[x y z elements is2d isplanet]=processmesh(md,[],options);
[elementnumbers datatype]=processdata(md,[1:md.mesh.numberofelements]',options);
edgecolor = getfieldvalue(options,'EdgeColor','black');

%plot
pos=getfieldvalue(options,'highlight',[]);
if is2d
	%plot mesh 
	A=elements(:,1); B=elements(:,2); C=elements(:,3);
	patch( 'Faces', [A B C], 'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor',edgecolor);

	%Highlight
	pos=getfieldvalue(options,'highlight',[]);
	A=elements(pos,1); B=elements(pos,2); C=elements(pos,3);
	patch( 'Faces', [A B C], 'Vertices', [x y z],'FaceVertexCData', [0.9 0.5 0.5],'FaceColor','flat','EdgeColor',edgecolor);
else
	if size(elements,2)==6, %prisms
		%plot mesh 
		A=elements(:,1); B=elements(:,2); C=elements(:,3); D=elements(:,4); E=elements(:,5); F=elements(:,6);
		patch( 'Faces', [A B C],  'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor',edgecolor);
		patch( 'Faces', [D E F],  'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor',edgecolor);
		patch( 'Faces', [A B E D],'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor',edgecolor);
		patch( 'Faces', [B E F C],'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor',edgecolor);
		patch( 'Faces', [C A D F],'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor',edgecolor);

		%Highlight
		A=elements(pos,1); B=elements(pos,2); C=elements(pos,3); D=elements(pos,4); E=elements(pos,5); F=elements(pos,6);
		patch( 'Faces', [A B C],  'Vertices', [x y z],'FaceVertexCData', [0.9 0.5 0.5],'FaceColor','flat','EdgeColor',edgecolor);
		patch( 'Faces', [D E F],  'Vertices', [x y z],'FaceVertexCData', [0.9 0.5 0.5],'FaceColor','flat','EdgeColor',edgecolor);
		patch( 'Faces', [A B E D],'Vertices', [x y z],'FaceVertexCData', [0.9 0.5 0.5],'FaceColor','flat','EdgeColor',edgecolor);
		patch( 'Faces', [B E F C],'Vertices', [x y z],'FaceVertexCData', [0.9 0.5 0.5],'FaceColor','flat','EdgeColor',edgecolor);
		patch( 'Faces', [C A D F],'Vertices', [x y z],'FaceVertexCData', [0.9 0.5 0.5],'FaceColor','flat','EdgeColor',edgecolor);
	elseif size(elements,2)==4, %tetras
		A=elements(:,1); B=elements(:,2); C=elements(:,3); D=elements(:,4);
		patch( 'Faces',[A B C],'Vertices', [x y z],'FaceVertexCData',zeros(size(x)),'FaceColor','none','EdgeColor',edgecolor);
		patch( 'Faces',[A B D],'Vertices', [x y z],'FaceVertexCData',zeros(size(x)),'FaceColor','none','EdgeColor',edgecolor);
		patch( 'Faces',[B C D],'Vertices', [x y z],'FaceVertexCData',zeros(size(x)),'FaceColor','none','EdgeColor',edgecolor);
		patch( 'Faces',[C A D],'Vertices', [x y z],'FaceVertexCData',zeros(size(x)),'FaceColor','none','EdgeColor',edgecolor);
		%Highlight
		A=elements(pos,1); B=elements(pos,2); C=elements(pos,3); D=elements(pos,4);
		patch( 'Faces',[A B C],'Vertices', [x y z],'FaceVertexCData',[0.9 0.5 0.5],'FaceColor','flat','EdgeColor',edgecolor);
		patch( 'Faces',[A B D],'Vertices', [x y z],'FaceVertexCData',[0.9 0.5 0.5],'FaceColor','flat','EdgeColor',edgecolor);
		patch( 'Faces',[B C D],'Vertices', [x y z],'FaceVertexCData',[0.9 0.5 0.5],'FaceColor','flat','EdgeColor',edgecolor);
		patch( 'Faces',[C A D],'Vertices', [x y z],'FaceVertexCData',[0.9 0.5 0.5],'FaceColor','flat','EdgeColor',edgecolor);
	else
		error('Not supported');
	end
end

%apply options
if ~exist(options,'highlight')
	disp('highlightelements warning : highlight option empty, not element highlighted');
end
options=addfielddefault(options,'title','Highlighted Elements');
options=addfielddefault(options,'colorbar',0);
applyoptions(md,[],options);
