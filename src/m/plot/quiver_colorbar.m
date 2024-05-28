function quiver_colorbar(quivers,options)
%QUIVER_COLORBAR - colorbar for quiver plots
%
%   Usage:
%      quiver_colorbar(quivers,options)

if  strcmpi(getfieldvalue(options,'colorbar','on'),'on'),

	%Create colorbar
	hcb=colorbar('peer',gca,'location','EastOutside');
	caxis([1 quivers.numcolors+1]);

	%build ticks
	ticklabel=cell(1,length(quivers.levels));
	for i=1:length(quivers.levels),
		ticklabel{i}=num2str(round_ice(quivers.levels(i),3));
	end
	tickpos=1:quivers.numcolors+1;

	%remove ticks if to many have been created
	proportion=round(length(quivers.levels)/4);
	if proportion>1,
		ticklabel=ticklabel(1:proportion:end);
		tickpos=tickpos(1:proportion:end);
	end

	%draw colorbar
	set(hcb,'YTickLabel',ticklabel,'YTick',tickpos);

	%position
	if exist(options,'colorbarpos'),
		set(hcb,'Position',getfieldvalue(options,'colorbarpos'));
	end
	%fontsize
	fontsize=getfieldvalue(options,'fontsize',14);
	set(hcb,'FontSize',fontsize);

	if exist(options,'colorbartitle'),
		backup=gca;
		axes(hcb);lab=title(getfieldvalue(options,'colorbartitle'));
		set(lab,'Rotation',getfieldvalue(options,'colorbartitlerotation',0));
		set(lab,'VerticalAlignment','bottom');
		axes(backup);
	end
end
