function areas=GetAreas(index,x,y,varargin)
%GETAREAS - compute areas or volumes of elements
%
%   compute areas of triangular elements or volumes 
%   of pentahedrons
%
%   Usage:
%      areas  =GetAreas(index,x,y);
%      volumes=GetAreas(index,x,y,z);
%
%   Examples:
%      areas  =GetAreas(md.mesh.elements,md.mesh.x,md.mesh.y);
%      volumes=GetAreas(md.mesh.elements,md.mesh.x,md.mesh.y,md.z);

%get number of elements and number of nodes
nels=size(index,1);
nods=length(x);
if nargin==4, z=varargin{1}; end

%some checks
if nargout~=1 | (nargin~=3 & nargin~=4),
	help GetAreas
	error('GetAreas error message: bad usage')
end
if ((length(y)~=nods) | (nargin==4 & length(z)~=nods)),
	error('GetAreas error message: x,y and z do not have the same length')
end
if max(index(:))>nods,
	error(['GetAreas error message: index should not have values above ' num2str(nods) ])
end
if (nargin==3 & size(index,2)~=3),
	error('GetAreas error message: index should have 3 columns for 2d meshes.')
end
if (nargin==4 & size(index,2)~=6),
	error('GetAreas error message: index should have 6 columns for 3d meshes.')
end

%initialization
areas=zeros(nels,1);
x1=x(index(:,1)); x2=x(index(:,2)); x3=x(index(:,3));
y1=y(index(:,1)); y2=y(index(:,2)); y3=y(index(:,3));

%compute the volume of each element
if nargin==3,
	%compute the surface of the triangle
	areas=(0.5*((x2-x1).*(y3-y1)-(y2-y1).*(x3-x1)));
else
	%V=area(triangle)*1/3(z1+z2+z3)
	thickness=mean(z(index(:,4:6)),2)-mean(z(index(:,1:3)),2);
	areas=(0.5*((x2-x1).*(y3-y1)-(y2-y1).*(x3-x1))).*thickness;
end
