function [vxout vyout] = interpJoughinCompositeGreenland(X,Y,ncpath),

%   - optional input argument: path to dataset IanGreenVel.mat

if nargin<3
	%data=load(['/u/astrid-r1b/morlighe/issmjpl/proj-morlighem/DatasetGreenland/Data/VelJoughin/IanGreenVel.mat']);
	filename = '/totten_1/ModelData/Greenland/VelJoughin/IanGreenVel.mat';
else
	filename = [ncpath '/IanGreenVel.mat']; 
end

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
x = x_m(id1x:id2x);
y = y_m(id1y:id2y);

vxout = InterpFromGrid(x,y,double(vx),X,Y);
vyout = InterpFromGrid(x,y,double(vy),X,Y);

if nargout==1,
	vxout = sqrt(vxout.^2+vyout.^2);
end
