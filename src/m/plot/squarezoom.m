function squarezoom()
%SQUAREZOOM - zoom on a part of a figure
%
%   Usage:
%      squarezoom()

disp('Click twice to define a square where you want to zoom. First click for upper left corner, second for lower right corner');
[x,y]=ginput(2);
dx=x(2)-x(1);
dy=y(1)-y(2);

if dx>dy,
	delta=dx-dy;
	xlim([x(1) x(2)]);
	ylim([y(2)-delta/2 y(1)+delta/2]);
else
	delta=dy-dx;
	xlim([x(1)-delta/2 x(2)+delta/2]);
	ylim([y(2) y(1)]);
end
