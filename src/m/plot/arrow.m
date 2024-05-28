function arrow(x0,y0,x1,y1,varargin)
%ARROW - plot arrow, using (x0,y0) and (x1,y1) as initial and end points. options can be specified.
%
%   Usage:
%      arrow(x1,y1,x2,y2,options)
%      where options is a lit of paired arguments of string OR enums
%      options can be: 
%            'ratio': default .5 (ratio headarrow/length)
%            'widthratio': default is 1/10 of length

%recover options
options=pairoptions(varargin{:});
ratio=getfieldvalue(options,'ratio',.5);
widthratio=getfieldvalue(options,'widthratio',.1);
color=getfieldvalue(options,'color','k');

%compute some values out of (x1,y1) and (x2,y2)
length=sqrt((x1-x0)^2+(y1-y0)^2);
width=length*widthratio;

%Build the two points Ap and Bp
x=zeros(2,1);
y=zeros(2,1);
x(1)=x0; y(1)=y0;
x(2)=x1; y(2)=y1;

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
p1=patch([A(1) B(1) C(1) D(1)],[A(2) B(2) C(2) D(2)],color);
set(p1,'EdgeColor',color); set(p1,'FaceColor',color);
p2=patch([E(1) F(1) G(1) H(1)],[E(2) F(2) G(2) H(2)],color);
set(p2,'EdgeColor',color); set(p2,'FaceColor',color);
