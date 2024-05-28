function areas=GetAreas3DTria(index,x,y,z,varargin)
%GETAREAS3DTRIA - compute areas of triangles with 3D coordinates 
%
%   Compute areas of triangles with 3D coordinates 
%
%   Usage:
%      areas=GetAreas3DTria(index,x,y,z);
%
%   Examples:
%      areas=GetAreas3DTria(md.mesh.elements,md.mesh.x,md.mesh.y,md.mesh.z);

%get number of elements and number of nodes
nels=size(index,1);
nods=length(x);

%some checks
if nargout~=1 | (nargin~=3 & nargin~=4),
	help GetAreas3DTria
	error('GetAreas3DTria error message: bad usage')
end
if ((length(y)~=nods) | (nargin==4 & length(z)~=nods)),
	error('GetAreas3DTria error message: x, y, and z do not have the same length')
end

if max(index(:))>nods,
	error(['GetAreas3DTria error message: index should not have values above ' num2str(nods) ])
end
if (nargin==4 & size(index,2)~=3),
	error('GetAreas3DTria error message: index should have 3 columns for 2d meshes')
end

%initialization
areas=zeros(nels,1);
x1=x(index(:,1)); x2=x(index(:,2)); x3=x(index(:,3));
y1=y(index(:,1)); y2=y(index(:,2)); y3=y(index(:,3));
z1=z(index(:,1)); z2=z(index(:,2)); z3=z(index(:,3));

%compute the volume of each element
if nargin==4,
	% area of triangles with 3D coordinates
	for i=1:nels
		m1=[x1(i) x2(i) x3(i); y1(i) y2(i) y3(i); 1 1 1];
		m2=[y1(i) y2(i) y3(i); z1(i) z2(i) z3(i); 1 1 1];
		m3=[z1(i) z2(i) z3(i); x1(i) x2(i) x3(i); 1 1 1];
		areas(i)=sqrt(det(m1)^2 + det(m2)^2 + det(m3)^2)/2;
	end
end
