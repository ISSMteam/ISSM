function expbox(filename)
%EXPBOX - Create an ARGUS file using two clicks
%
%   Two clicks on a plot are used to generate a rectangular box
%   This box is written in EXP format on filename
%
%   Usage:
%      expbox(filename)

%check
if exist(filename,'file'),
	choice=input(['A file ' filename ' already exists, do you want to modify it? (y/n)'],'s');
	if ~strcmpi(choice,'y'),
		disp('no modification done ... exiting');
		return
	end
end

%Get points
disp('Click twice to define a rectangular domain. First click for upper left corner, second for lower right corner');
[x,y]=ginput(2);

x1=x(1);
x2=x(2);
x3=x2;
x4=x1;

y1=y(1);
y2=y1;
y3=y(2);
y4=y3;

%Build Exp structure
A=struct();
A.nods=5;
A.density=1;
A.x=[x1 x2 x3 x4 x1]';
A.y=[y1 y2 y3 y4 y1]';

%Write structure
expwrite(A,filename);
