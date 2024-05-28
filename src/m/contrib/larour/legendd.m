function legendd(varargin)

	options=pairoptions(varargin{:});

	%retrieve arguments: 
	x=getfieldvalue(options,'x',0); 
	y=getfieldvalue(options,'y',0);
	w=getfieldvalue(options,'w',1);
	h=getfieldvalue(options,'h',1);
	facecolor=getfieldvalue(options,'FaceColor','w');
	edgecolor=getfieldvalue(options,'EdgeColor','k');
	strings=getfieldvalue(options,'strings',{});
	colors=getfieldvalue(options,'colors',{});
	fontsize=getfieldvalue(options,'FontSize',12);
	linewidth=getfieldvalue(options,'LineWidth',2);

	hold on;
	rectangle('Position',[x,y,w,h],'FaceColor',facecolor,'EdgeColor',edgecolor);
	
	nl=length(strings);
	for i=1:nl,
		l=line([x+w/6 x+w/3],[y+(nl+1-i)*h/(nl+1) y+(nl+1-i)*h/(nl+1)]);
		set(l,'Color',colors{i});
		set(l,'LineWidth',linewidth);
		text(x+1.3*w/3,y+(nl+1-i)*h/(nl+1),strings{i},'FontSize',fontsize);
	end

