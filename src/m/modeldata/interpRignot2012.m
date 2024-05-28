function [vxout vyout]= interpRignot2012(X,Y),

filename = '/totten_1/ModelData/Greenland/VelMouginot/RignotGreenland2012Vel.mat';


%Figure out what subset of the matrix should be read
load(filename,'x','y');
velfile = matfile(filename);

offset=2;

xmin=min(X(:)); xmax=max(X(:));
posx=find(x<=xmax);
id1x=max(1,find(x>=xmin,1)-offset);
id2x=min(numel(x),posx(end)+offset);

ymin=min(Y(:)); ymax=max(Y(:));
%posy=find(y>=ymin);
%id1y=max(1,find(y<=ymax,1)-offset);
%id2y=min(numel(y),posy(end)+offset);
posy=find(y<=ymax);
id1y=max(1,find(y>=ymin,1)-offset);
id2y=min(numel(y),posy(end)+offset);

vx = velfile.vx(id1y:id2y,id1x:id2x);
vy = velfile.vy(id1y:id2y,id1x:id2x);
x = x(id1x:id2x);
y = y(id1y:id2y);

%load(filename);
vxout = InterpFromGrid(x,y,double(vx),X,Y);
vyout = InterpFromGrid(x,y,double(vy),X,Y);

if nargout==1,
	vxout = sqrt(vxout.^2+vyout.^2);
end
