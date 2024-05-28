function [vxout vyout]= interpMouginotAnt2016(X,Y),

%read data
switch (oshostname()),
	case {'ronne'}
		filename = '/home/ModelData/Antarctica/MouginotVel/vel_ant_5Apr2016.mat';
	case {'thwaites','murdo','astrid'}
		filename = '/u/astrid-r1b/ModelData/RignotAntarcticaVelMosaic450m/vel_ant_5Apr2016.mat';
	otherwise
		error('hostname not supported yet');
end

%Figure out what subset of the matrix should be read
load(filename,'x','y');
velfile = matfile(filename);

offset=2;

xmin=min(X(:)); xmax=max(X(:));
posx=find(x<=xmax);
id1x=max(1,find(x>=xmin,1)-offset);
id2x=min(numel(x),posx(end)+offset);

if y(2)-y(1)<0
	ymin=min(Y(:)); ymax=max(Y(:));
	posy=find(y>=ymin);
	id1y=max(1,find(y<=ymax,1)-offset);
	id2y=min(numel(y),posy(end)+offset);
else
	ymin=min(X(:)); ymax=max(X(:));
	posy=find(y<=ymax);
	id1y=max(1,find(y>=ymin,1)-offset);
	id2y=min(numel(y),posy(end)+offset);
end

vx = velfile.vx(id1y:id2y,id1x:id2x);
vy = velfile.vy(id1y:id2y,id1x:id2x);
x = x(id1x:id2x);
y = y(id1y:id2y);

vxout = InterpFromGrid(x,y,double(vx),X,Y);
vyout = InterpFromGrid(x,y,double(vy),X,Y);

if nargout==1,
	vxout = sqrt(vxout.^2+vyout.^2);
end
