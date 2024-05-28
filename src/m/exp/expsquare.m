function expsquare(filename)
%EXPBOX - Create a ARGUS file using to clicks
%
%   Two clicks on a plot are used to generate a square box
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
disp('Click twice to define a square domain. First click for upper left corner, second for lower right corner');
[x,y]=ginput(2);

xmiddle=mean(x);
ymiddle=mean(y);

x1=x(1); y1=y(1);
x3=x(2); y3=y(2);

Diag=[x1-xmiddle;y1-ymiddle];

Vector=[xmiddle;ymiddle]+[-Diag(2);Diag(1)];
x2=Vector(1); y2=Vector(2);

Vector=[xmiddle;ymiddle]-[-Diag(2);Diag(1)];
x4=Vector(1); y4=Vector(2);

%Build Exp structure
A=struct();
A.nods=5;
A.density=1;
A.x=[x1 x2 x3 x4 x1]';
A.y=[y1 y2 y3 y4 y1]';

%Write structure
expwrite(A,filename);
