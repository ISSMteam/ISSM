function [vxout vyout] = interpJoughinMosaic(X,Y),

switch oshostname(),
	case {'ronne'}
		filename = '/home/ModelData/Greenland/VelJoughin/IanGreenVel.mat';
	case {'totten'}
		filename = '/totten_1/ModelData/Greenland/VelJoughin/IanGreenVel.mat';
	otherwise
		error('machine not supported yet');
end
verbose = 1;

%Figure out what subset of the matrix should be read
load(filename,'x_m','y_m');
velfile = matfile(filename);

offset=2;

xmin=min(X(:)); xmax=max(X(:));
posx=find(x_m<=xmax);
id1x=max(1,find(x_m>=xmin,1)-offset);
id2x=min(numel(x_m),posx(end)+offset);

ymin=min(Y(:)); ymax=max(Y(:));
posy=find(y_m>=ymin);
id1y=max(1,find(y_m<=ymax,1)-offset);
id2y=min(numel(y_m),posy(end)+offset);

vx = velfile.vx(id1y:id2y,id1x:id2x);
vy = velfile.vy(id1y:id2y,id1x:id2x);
x_m = x_m(id1x:id2x);
y_m = y_m(id1y:id2y);

%load(filename);
vxout = InterpFromGrid(x_m,y_m,vx,X,Y);
vyout = InterpFromGrid(x_m,y_m,vy,X,Y);
