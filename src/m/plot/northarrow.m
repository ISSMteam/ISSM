function northarrow(structure)
%NORTHARROW - overlay an arrow pointing north on the current plot
%
%   Usage:
%      northarrow(structure)

%Go through structure and fill missing arguments
if length(structure)<3
	error('plotmodel error message: the position or the length of the North arrow is missing');
elseif length(structure)==3
	structure(4)=0.5; %default ratio headarrow/length
	structure(5)=structure(3)/10; %default width =length/10
elseif length(structure)==4
	structure(5)=structure(3)/10; %default width =length/10
elseif length(structure)==5
	structure(6)=16; %default fontsize
elseif length(structure)>6
	error('plotmodel error message: to many input arguments for northarrow: [x0 y0 length [ratio width fontsize]]');
end

%retrieve north arrow parameters
x0=structure(1);
y0=structure(2);
lengtharrow=structure(3);
ratio=structure(4);
width=structure(5);
fontsize=structure(6);

%Figure out angle to point towards north
ang=atan2(y0,x0);

%Build the two points Ap and Bp
x=zeros(2,1);
y=zeros(2,1);
x(1)=x0;
y(1)=y0;

x(2)=x(1)+lengtharrow*cos(ang);
y(2)=y(1)+lengtharrow*sin(ang);

Ap=[x(1)
   y(1)];
Bp=[x(2)
   y(2)];

%Build arrowhead first
ang2=150*2*pi/360;
rotation=[cos(ang2), sin(ang2); -sin(ang2), cos(ang2)];

E=ratio*rotation*(Bp-Ap)+Bp;
F=Bp;
G=ratio*rotation'*(Bp-Ap)+Bp;
H=Bp/4+E*3/8+G*3/8;

%Build rectangle
u=Bp-Ap;
alpha=atan2(u(2),u(1));

A=Ap-[-width/2*sin(alpha)
   width/2*cos(alpha)];
 B=H-[-width/2*sin(alpha)
   width/2*cos(alpha)];
C=H+[-width/2*sin(alpha)
   width/2*cos(alpha)];
D=Ap+[-width/2*sin(alpha)
   width/2*cos(alpha)];

%Plot arrow
hold on
p1=patch([A(1) B(1) C(1) D(1)],[A(2) B(2) C(2) D(2)],'Black');
p2=patch([E(1) F(1) G(1) H(1)],[E(2) F(2) G(2) H(2)],'Black');

%Text North
xN=max([A(1) D(1) E(1) F(1) G(1)])+ratio/3*abs(lengtharrow);
yN=mean([A(2) F(2) H(2)]);
text(xN,yN,'North','FontSize',fontsize,'FontWeight','b');
