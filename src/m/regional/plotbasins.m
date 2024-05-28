%display all the domain outlines in a directory

basins=listfiles;

hold on
for i=1:length(basins), 
	%check whether this is a .exp file
	basin=basins{i};
	if strcmpi(basin(end-3:end),'.exp'),

		contour=expread(basin,0);
		x=contour(1).x;
		y=contour(1).y;
		x0=mean(x); y0=mean(y);
		text(x0,y0,basin(1:end-4),'Fontsize',14);
		expdisp(basin);
	end
end
