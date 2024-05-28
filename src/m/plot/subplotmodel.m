function ha=subplotmodel(nlines,ncols,num,options)
%SUBPLOTMODEL -  tight subplot that includes margins
%
%   Usage:
%      h=subplotmodel(nlines,ncols,i,options)

box=getfieldvalue(options,'box','off');
if strcmpi(box,'on'),
	visible='on';
else
	visible='off';
end

%Do we specify axes?
if exist(options,'axes')
	ax = getfieldvalue(options,'axes');
	ha = axes('Units','normalized', ...
		'Position',ax,'XTickLabel','','YTickLabel','',...
		'Visible',visible,'box',box);
	return
end

%Regular subplot
if ~exist(options,'tightsubplot')
	if exist(options,'asymsubplot')
		id=getfieldvalue(options,'asymsubplot',num);
		subplot(nlines,ncols,id);
		return
	else
		subplot(nlines,ncols,num);
		return;
	end
end

gap     = getfieldvalue(options,'gap',[.01 .01]);
hmargin = getfieldvalue(options,'hmargin',[.01 .01]);
vmargin = getfieldvalue(options,'vmargin',[.01 .01]);

height = (1-sum(vmargin)-(nlines-1)*gap(1))/nlines; 
width  = (1-sum(hmargin)-(ncols-1)*gap(2))/ncols;
ymin   = 1-vmargin(2)-height; 

for i = 1:nlines
	xmin = hmargin(1);
	for j = 1:ncols
		if(((i-1)*ncols+j)==num)
			ha = axes('Units','normalized', ...
				'Position',[xmin ymin width height],'XTickLabel','','YTickLabel','',...
				'Visible',visible,'box',box);
			return
		end
		xmin = xmin+width+gap(2);
	end
	ymin = ymin-height-gap(1);
end

%Activate new axes
axes(ha);
