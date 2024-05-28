function h2 = plot_icefront(md,options,width,i,data)
%PLOT_ICEFRONT - plot segment on neumann BC
%
%   Usage:
%      plot_icefront(md,options,width,i);
%
%   See also: PLOTMODEL

%plot mesh boundaries
subplot(width,width,i); 

%process mesh and data
[x y z elements is2d isplanet]=processmesh(md,[],options);
ice=(md.mask.ice_levelset<0);
noice=(md.mask.ice_levelset>=0);
zeroice=(md.mask.ice_levelset==0);
elementice=sum(ice(md.mesh.elements),2);
elementnoice=sum(noice(md.mesh.elements),2);
elementzeroice=sum(zeroice(md.mesh.elements),2);

if is2d,
	icefront=(elementice & elementnoice) & ~(elementice==2 & elementzeroice);

	%plot mesh
	A=elements(:,1); B=elements(:,2); C=elements(:,3); 
	h1=patch( 'Faces', [A B C], 'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
	hold on;

	%highlight elements on neumann
	pos=find(icefront);
	A=elements(pos,1); B=elements(pos,2); C=elements(pos,3); 
	h2=patch( 'Faces', [A B C], 'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','blue','EdgeColor','black');
	hold on;

	%Plot zero ice_levelset line

else
	icefront=(elementice & elementnoice) & ~(elementice==4 & elementzeroice);

	%plot mesh
	A=elements(:,1); B=elements(:,2); C=elements(:,3); D=elements(:,4); E=elements(:,5); F=elements(:,6);
	h1=patch( 'Faces', [A B C],'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
	patch( 'Faces', [D E F],  'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
	patch( 'Faces', [A B E D],'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
	patch( 'Faces', [B E F C],'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
	patch( 'Faces', [C A D F],'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
	hold on;

	%highlight elements on neumann
	pos=find(icefront);
	A=elements(pos,1); B=elements(pos,2); C=elements(pos,3); D=elements(pos,4); E=elements(pos,5); F=elements(pos,6);
	h2=patch( 'Faces', [A B C],'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','blue','EdgeColor','black');
	patch( 'Faces', [D E F],  'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','blue','EdgeColor','black');
	patch( 'Faces', [A B E D],'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','blue','EdgeColor','black');
	patch( 'Faces', [B E F C],'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','blue','EdgeColor','black');
	patch( 'Faces', [C A D F],'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','blue','EdgeColor','black');
end

%legend (disable warnings)
warning off
legend([h2],'element on ice front')
warning on

%apply options
options=addfielddefault(options,'title','Neumann boundary conditions');
options=addfielddefault(options,'colorbar',0);
applyoptions(md,[],options);
