function plot_elementstype(md,options,width,i)
%PLOT_ELEMENTSTYPE - plot elements type
%
%   Usage:
%      plot_elementstype(md,options,width,i);
%
%   See also: PLOTMODEL

%process data and model
[x y z elements is2d isplanet]=processmesh(md,[],options);
[data datatype]=processdata(md,md.flowequation.element_equation,options);

%edgecolor?
edgecolor=getfieldvalue(options,'edgecolor','none');

%plot
subplot(width,width,i);
p = [];

if is2d
	for i=0:9
		pos=find(data==i);
		A=elements(pos,1); B=elements(pos,2); C=elements(pos,3); 
		pnew=patch( 'Faces', [A B C], 'Vertices', [x y z],'CData',i,'FaceColor','flat','EdgeColor',edgecolor);
		p = [p;pnew];
	end
else
	for i=0:9
		pos=find(data==i);
		A=elements(pos,1); B=elements(pos,2); C=elements(pos,3); D=elements(pos,4); E=elements(pos,5); F=elements(pos,6);
		pnew = patch( 'Faces', [A B C],  'Vertices', [x y z],'CData', i,'FaceColor','flat','EdgeColor',edgecolor);
		p = [p;pnew];
		patch( 'Faces', [D E F],  'Vertices', [x y z],'CData', i,'FaceColor','flat','EdgeColor',edgecolor);
		patch( 'Faces', [A B E D],'Vertices', [x y z],'CData', i,'FaceColor','flat','EdgeColor',edgecolor);
		patch( 'Faces', [B E F C],'Vertices', [x y z],'CData', i,'FaceColor','flat','EdgeColor',edgecolor);
		patch( 'Faces', [C A D F],'Vertices', [x y z],'CData', i,'FaceColor','flat','EdgeColor',edgecolor);
	end
end
legend(p,'None','SIA','SSA','L1L2','MOLHO','HO',...
		'SSAHO','FS','SSAFS','HOFS');

%apply options
options=addfielddefault(options,'title','Elements type');
options=addfielddefault(options,'colorbar',0);
applyoptions(md,[],options);
