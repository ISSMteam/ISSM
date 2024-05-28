function md=squaremesh(md,Lx,Ly,nx,ny,computeconnectivity,orientation)
%SQUAREMESH - create a structured square mesh 
%
%   This script will generate a structured square mesh
%   Lx and Ly are the dimension of the domain (in meters)
%   nx anx ny are the number of nodes in the x and y direction
%   The coordinates x and y returned are in meters.
%
%   Usage:
%      [md]=squaremesh(md,Lx,Ly,nx,ny)

%process options
if nargin == 5
	computeconnectivity = 1;
	orientation = 1;
elseif nargin==6
	orientation = 1;
end

%get number of elements and number of nodes
nel=(nx-1)*(ny-1)*2;
nods=nx*ny;

%prepare coordinates of vertices
x = repmat(linspace(0,Lx,nx),[ny 1]);
x = reshape(x,[nx*ny 1]);
y = repmat(linspace(0,Ly,ny)',[1 nx]);
y = reshape(y,[nx*ny 1]);

%do first column of elements first
nels1 = 2*(ny-1);
index = ones(nels1,3);
%First column
index(2:2:nels1,1) = 2:ny;
index(3:2:nels1,1) = 2:ny-1;
%2d column
index(1:2:nels1,2) = ny+1:2*ny-1;
index(2:2:nels1,2) = ny+1:2*ny-1;
%3rd column
index(1:2:nels1,3) = 2:ny;
index(2:2:nels1,3) = ny+2:2*ny;

%Now copy column and offset with ny, nx times
index = repmat(index,[nx-1 1]);
offset = repmat([0:ny:(nx-2)*ny],[nels1 1]);
offset = reshape(offset,[(nx-1)*nels1,1]);
offset = repmat(offset,[1,3]);
index = index + offset;

%create segments
segments=zeros(2*(nx-1)+2*(ny-1),3);
%left edge:
segments(1:ny-1,:)=[[2:ny]' [1:ny-1]' 2*[1:ny-1]'-1];
%right edge:
segments(ny:2*(ny-1),:)=[[ny*(nx-1)+1:nx*ny-1]' [ny*(nx-1)+2:nx*ny]' 2*[(ny-1)*(nx-2)+1:(nx-1)*(ny-1)]'];
%front edge:
segments(2*(ny-1)+1:2*(ny-1)+(nx-1),:)=[[2*ny:ny:ny*nx]' [ny:ny:ny*(nx-1)]' [2*(ny-1):2*(ny-1):2*(nx-1)*(ny-1)]'];
%back edge
segments(2*(ny-1)+(nx-1)+1:2*(nx-1)+2*(ny-1),:)=[[1:ny:(nx-2)*ny+1]' [ny+1:ny:ny*(nx-1)+1]' [1:2*(ny-1):2*(nx-2)*(ny-1)+1]'];

%Do we need to change the orientation of the diagonals?
if orientation==1
	%nothing to do
elseif orientation==2
	%switch diagonals
	indexold = index;
	index(1:2:end,3) = indexold(2:2:end,3);
	index(2:2:end,2) = indexold(1:2:end,1);
elseif orientation==3
	%alternate diagonals
	indexold = index;
	index(1:4:end,3) = indexold(2:4:end,3);
	index(2:4:end,2) = indexold(1:4:end,1);
else
	error('not supported');
end

%plug coordinates and nodes
md.mesh=mesh2d();
md.mesh.x=x;
md.mesh.y=y;
md.mesh.numberofvertices=nods;
md.mesh.vertexonboundary=zeros(nods,1);md.mesh.vertexonboundary(segments(:,1:2))=1;

%plug elements
md.mesh.elements=index;
md.mesh.segments=segments;
md.mesh.numberofelements=nel;

%Now, build the connectivity tables for this mesh.
if computeconnectivity
	md.mesh.vertexconnectivity=NodeConnectivity(md.mesh.elements,md.mesh.numberofvertices);
	md.mesh.elementconnectivity=ElementConnectivity(md.mesh.elements,md.mesh.vertexconnectivity);
end
