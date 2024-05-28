function plot_unit(x,y,z,elements,data,is2d,isplanet,datatype,options)
%PLOT_UNIT - unit plot, display data
%
%   Usage:
%      plot_unit(x,y,z,elements,data,is2d,isplanet,datatype,options);
%
%   See also: PLOTMODEL, PLOT_MANAGER

%edgecolor
edgecolor=getfieldvalue(options,'edgecolor','none');

switch datatype,

	%element plot
	case 1,

		pos=find(~isnan(data)); %needed for element on water
		if size(elements,2)==6, %prisms
			A=elements(pos,1); B=elements(pos,2); C=elements(pos,3); D=elements(pos,4); E=elements(pos,5); F=elements(pos,6);
			patch( 'Faces', [A B C],  'Vertices', [x y z],'CData', data(pos),'FaceColor','flat','EdgeColor',edgecolor);
			patch( 'Faces', [D E F],  'Vertices', [x y z],'CData', data(pos),'FaceColor','flat','EdgeColor',edgecolor);
			patch( 'Faces', [A B E D],'Vertices', [x y z],'CData', data(pos),'FaceColor','flat','EdgeColor',edgecolor);
			patch( 'Faces', [B E F C],'Vertices', [x y z],'CData', data(pos),'FaceColor','flat','EdgeColor',edgecolor);
			patch( 'Faces', [C A D F],'Vertices', [x y z],'CData', data(pos),'FaceColor','flat','EdgeColor',edgecolor);
		elseif size(elements,2)==4, %tetras
			A=elements(:,1); B=elements(:,2); C=elements(:,3); D=elements(:,4);
			patch( 'Faces',[A B C],'Vertices', [x y z],'CData',data(pos),'FaceColor','flat','EdgeColor',edgecolor);
			patch( 'Faces',[A B D],'Vertices', [x y z],'CData',data(pos),'FaceColor','flat','EdgeColor',edgecolor);
			patch( 'Faces',[B C D],'Vertices', [x y z],'CData',data(pos),'FaceColor','flat','EdgeColor',edgecolor);
			patch( 'Faces',[C A D],'Vertices', [x y z],'CData',data(pos),'FaceColor','flat','EdgeColor',edgecolor);
		else
			A=elements(pos,1); B=elements(pos,2); C=elements(pos,3); data=data(pos);
			patch( 'Faces', [A B C], 'Vertices', [x y z],'CData', data,'FaceColor','flat','EdgeColor',edgecolor);
			
			%mask value NaN, plot white faces.
			if getfieldvalue(options,'maskwhite',0),
				pos2=find(data==getfieldvalue(options,'maskvalue',NaN));
				if ~isempty(pos2),
					A=elements(pos(pos2),1); B=elements(pos(pos2),2); C=elements(pos(pos2),3); data=data(pos2);
					patch( 'Faces', [A B C], 'Vertices', [x y z],'FaceColor','w','EdgeColor',edgecolor);
				end
			end
		end

		if is2d,
		end

	%node plot
	case 2,

		if size(elements,2)==6, %prisms
			A=elements(:,1); B=elements(:,2); C=elements(:,3); D=elements(:,4); E=elements(:,5); F=elements(:,6);
			patch( 'Faces', [A B C],  'Vertices', [x y z],'FaceVertexCData', data(:),'FaceColor','interp','EdgeColor',edgecolor);
			patch( 'Faces', [D E F],  'Vertices', [x y z],'FaceVertexCData', data(:),'FaceColor','interp','EdgeColor',edgecolor);
			patch( 'Faces', [A B E D],'Vertices', [x y z],'FaceVertexCData', data(:),'FaceColor','interp','EdgeColor',edgecolor);
			patch( 'Faces', [B E F C],'Vertices', [x y z],'FaceVertexCData', data(:),'FaceColor','interp','EdgeColor',edgecolor);
			patch( 'Faces', [C A D F],'Vertices', [x y z],'FaceVertexCData', data(:),'FaceColor','interp','EdgeColor',edgecolor);
		elseif size(elements,2)==4, %tetras
			A=elements(:,1); B=elements(:,2); C=elements(:,3); D=elements(:,4);
			patch( 'Faces',[A B C],'Vertices', [x y z],'FaceVertexCData',data(:),'FaceColor','interp','EdgeColor',edgecolor);
			patch( 'Faces',[A B D],'Vertices', [x y z],'FaceVertexCData',data(:),'FaceColor','interp','EdgeColor',edgecolor);
			patch( 'Faces',[B C D],'Vertices', [x y z],'FaceVertexCData',data(:),'FaceColor','interp','EdgeColor',edgecolor);
			patch( 'Faces',[C A D],'Vertices', [x y z],'FaceVertexCData',data(:),'FaceColor','interp','EdgeColor',edgecolor);
		else
			A=elements(:,1); B=elements(:,2); C=elements(:,3); 
			patch( 'Faces', [A B C], 'Vertices', [x y z],'FaceVertexCData', data(:),'FaceColor','interp','EdgeColor',edgecolor);
		end

	%quiver plot
	case 3,

		if is2d,
			plot_quiver(x,y,data(:,1),data(:,2),options);
		else
			plot_quiver3(x,y,z,data(:,1),data(:,2),data(:,3),options);
		end

	%edge data
	case 6
		A=elements(:,1); B=elements(:,2); C=elements(:,3); 
		patch('Faces', [A B C],'Vertices', [x y z],'FaceVertexCData', data(:),'FaceColor','interp','EdgeColor',edgecolor);
	otherwise,
		error(['case ' num2str(datatype) ' not supported']);

	end
end
