function plot_quiver(x,y,u,v,options),
%PLOT_QUIVER - quiver plot with colors
%
%   to be perfected tomorrow
%
%   Usage:
%      plot_quiver(x,y,u,v,options)
%
%   Example:
%      plot_quiver(md.mesh.x,md.mesh.y,md.initialization.vx,md.initialization.vy,options);

%process fields
[quivers,palette]=quiver_process(x,y,u,v,options);

%loop over the number of colors
hold on
h=[];
for i=1:quivers.numcolors
	pos=find(quivers.colorind==i);
	if ~isempty(pos),
		hprime=quiver(quivers.x(pos),quivers.y(pos),quivers.u(pos),quivers.v(pos),...
			'Color',palette(i,:),'ShowArrowHead','on','AutoScale','off');
		h=[h;hprime];
	end
end

%take care of colorbar
quiver_colorbar(quivers,options);
