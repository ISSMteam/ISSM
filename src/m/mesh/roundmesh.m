function md=roundmesh(md,radius,resolution,varargin)
%ROUNDMESH - create an unstructured round mesh 
%
%   This script will generate an unstructured round mesh
%   - radius     : specifies the radius of the circle in meters
%   - resolution : specifies the resolution in meters
%
%   Usage:
%      md=roundmesh(md,radius,resolution)
%      md=roundmesh(md,radius,resolution,'domain.exp')

%First we have to create the domain outline 
if nargin>=4
	expname = varargin{1};
else
	expname = [tempname() '.exp'];
end

%Get number of points on the circle
pointsonedge=floor((2.*pi*radius) / resolution)+1; %+1 to close the outline

%Calculate the Cartesian coordinates of the points
theta=linspace(0,2*pi,pointsonedge)';
x_list=roundsigfig(radius*cos(theta),12);
y_list=roundsigfig(radius*sin(theta),12);
A=struct('x',x_list,'y',y_list,'density',1.);
expwrite(A,expname);

%Call mesher
md=triangle(md,expname,resolution);
%md=bamg(md,'domain','RoundDomainOutline.exp','hmin',resolution);

%move the closest node to the center
[minimum pos]=min(md.mesh.x.^2+md.mesh.y.^2);
md.mesh.x(pos)=0.;
md.mesh.y(pos)=0.;

%delete domain
if nargin<4
	delete(expname);
end
end

function x=roundsigfig(x,n)

digits=ceil(log10(abs(x)));
x=x./10.^digits;
x=round(x.*10.^n)./10.^n;
x=x.*10.^digits;

pos=find(isnan(x));
x(pos)=0.;

end
