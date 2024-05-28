function curvedarrow(centerx,centery,distance,angle,secondangle,varargin)
%CURVEDARROW - plot curved arrow, where curvature center is on (centerx,centery), radius is 'distance', 
%   and angle determines arc length (angle in degrees). secondangle determines by how much the arc 
%   generated is rotated along the [0,0,1] z axis.30,
%
%   Usage:
%      curevedarrow(x1,y1,r,30,30,options)
%      where options is a lit of paired arguments of string OR enums
%      options can be: 
%            'ratio': default .5 (ratio headarrow/length)
%            'widthratio': default is 1/10 of length
%            'width': if you want to specify an absolute width

	%recover options
	options=pairoptions(varargin{:});
	ratio=getfieldvalue(options,'ratio',.1);
	arrowlength=getfieldvalue(options,'arrowlength',1);
	color=getfieldvalue(options,'color','k');

	%transform angle in radians
	angle=angle/180*pi;
	nsteps=10;

	%compute some values out of (x1,y1) and (x2,y2)
	length=distance*angle;

	if exist(options,'widthratio'),
		widthratio=getfieldvalue(options,'widthratio');
		width=length*widthratio;
	else if exist(options,'width'),
		width=getfieldvalue(options,'width');
	else 
		widthratio=.1;
		width=length*widthratio;
	end

	%buidl the arrow itself: 
	A=[centerx+distance, centery];
	B=[centerx+distance+width, centery];
	BC=[centerx+cos(0:angle/nsteps:angle)*(distance+width); centery+sin(0:angle/nsteps:angle)*(distance+width)]';
	C=[centerx+cos(angle)*(distance+width), centery+sin(angle)*(distance+width)];
	D=[centerx+cos(angle)*(distance), centery+sin(angle)*(distance)];
	DA=[centerx+cos(angle:-angle/nsteps:0)*(distance); centery+sin(angle:-angle/nsteps:0)*(distance)]';

	%Plot arrow
	hold on
	p1=patch([A(1) B(1) BC(:,1)' C(1) D(1) DA(:,1)' ],[A(2) B(2) BC(:,2)' C(2) D(2) DA(:,2)'],color);
	set(p1,'EdgeColor',color); set(p1,'FaceColor',color);

	%Build arrowhead 
	E=D+2/3*(D-C);
	F=C+2/3*(C-D);

	n=(F-E)/norm(F-E,2);
	m=[-n(2) n(1)];
	if(angle<0)
		m=-m;
	end

	if exist(options,'arrowlength'),
		d=arrowlength;
	else
		d=abs((distance*angle)*ratio);
	end

	%G is d distance from middle of E and F: 
	G=(E+F)/2+d*m;

	p2=patch([E(1) F(1) G(1)],[E(2) F(2) G(2)],color);
	set(p2,'EdgeColor',color); set(p2,'FaceColor',color);

	%now rotate our arrow: 
	rotate(p1,[0 0  1],secondangle);
	rotate(p2,[0 0 1],secondangle);

end
