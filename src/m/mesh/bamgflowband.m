function md=bamgflowband(md,x,surf,base,varargin);
%BAMGFLOWBAND - create flowband mesh with bamg
%
%   Usage:
%      md=bamgflowband(md,x,surf,base,OPTIONS)
%
%      surf and bed are the surface elevation and base for each x provided
%      x must be increasing
%      OPTIONS are bamg options
%
%   Example:
%      x =[1:100:3000];
%      h=linspace(1000,300,numel(x));
%      b=-917/1023*h;
%      md=bamgflowband(model,b+h,b,'hmax',80,'vertical',1,'Markers',m);

%Write expfile with domain outline
A=struct();
A.x=[x;flipud(x);x(1)];
A.y=[base;flipud(surf);base(1)];
A.nods = numel(A.x);

%markers:
m                          = ones(numel(A.x)-1,1); % base        = 1
m(numel(x))                = 2;                    % right side  = 2
m(numel(x)+1:2*numel(x)-1) = 3;                    % top surface = 3
m(2*numel(x))              = 4;                    % left side   = 4

%mesh domain
md=bamg(model(),'domain',A,'vertical',1,'Markers',m,varargin{:});

%Deal with vertices on bed
md.mesh.vertexonbase=zeros(md.mesh.numberofvertices,1);
md.mesh.vertexonbase(find(vertexflags(md.mesh,1)))=1;
md.mesh.vertexonsurface=zeros(md.mesh.numberofvertices,1);
md.mesh.vertexonsurface(find(vertexflags(md.mesh,3)))=1;
