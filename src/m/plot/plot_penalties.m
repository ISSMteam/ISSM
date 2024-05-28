function plot_penalties(md,options,width,i)
%PLOT_PENALTIES - plot penalties
%
%   Usage:
%      plot_penalties(md,options,width,i);
%
%   See also: PLOTMODEL

%process data and model
[x y z elements is2d isplanet]=processmesh(md,[],options);

%plot mesh penalties
subplot(width,width,i); 

%units
if exist(options,'unit'),
	unit=getfieldvalue(options,'unit');
	x=x*unit;
	y=y*unit;
	z=z*unit;
end

if dimension(md.mesh)~=3,
	error('no penalties to plot for ''2d'' model');
elseif isempty(md.penalties),
	disp('no penalty applied in this model');
	return;
else
	%plot mesh
	A=elements(:,1); B=elements(:,2); C=elements(:,3); D=elements(:,4); E=elements(:,5); F=elements(:,6);
	patch( 'Faces', [A B C],  'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
	patch( 'Faces', [D E F],  'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
	patch( 'Faces', [A B E D],'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
	patch( 'Faces', [B E F C],'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
	patch( 'Faces', [C A D F],'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');

	hold on;
	for (i=1:size(md.penalties,1)),
		P1=plot3(x(md.penalties(i,1)),y(md.penalties(i,1)),z(md.penalties(i,1)),'ro','MarkerSize',15,'MarkerFaceColor','r');
		P2=plot3(x(md.penalties(i,:)),y(md.penalties(i,:)),z(md.penalties(i,:)),'bo-','LineWidth',2,'MarkerSize',8,'MarkerFaceColor','b');
	end
	legend([P1 P2],'SSA''s penalized nodes','HO''s penalized nodes');
end

%apply options
options=addfielddefault(options,'title','Penalties');
options=addfielddefault(options,'colorbar',0);
applyoptions(md,[],options);
