function md=rifttipsrefine(md,filename,resolution,circleradius)
%RIFTTIPSREFINE - refine mesh near rift tips
%
%   Usage:
%      md=rifttipsrefine(md,filename,resolution,circleradius);

numberofnodes=50;

%take rifts, and create refinement circles around tips
rifts=expread(filename);

!echo -n "" > Circles.exp
for i=1:length(rifts),
	tip1=[rifts(i).x(1) rifts(i).y(1)];
	tip2=[rifts(i).x(end) rifts(i).y(end)];
	%create circle around tip
	expcreatecircle('Circle1.exp',tip1(1),tip1(2),circleradius,numberofnodes);
	expcreatecircle('Circle2.exp',tip2(1),tip2(2),circleradius,numberofnodes);
	!cat Circles.exp Circle1.exp Circle2.exp > Circles2.exp
	!mv Circles2.exp Circles.exp
	!rm -rf Circle1.exp Circle2.exp
end

md=meshexprefine(md,'Circles.exp',resolution);

system('rm -rf Circles.exp');
