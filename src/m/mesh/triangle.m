function md=triangle(md,domainname,varargin)
%TRIANGLE - create model mesh using the triangle package
%
%   This routine creates a model mesh using Triangle and a domain outline, to within a certain resolution
%   where md is a @model object, domainname is the name of an Argus domain outline file, 
%   and resolution is a characteristic length for the mesh (same unit as the domain outline
%   unit). Riftname is an optional argument (Argus domain outline) describing rifts.
%
%   Usage:
%      md=triangle(md,domainname,resolution)
%   or md=triangle(md,domainname,riftname, resolution)
%
%   Examples:
%      md=triangle(md,'DomainOutline.exp',1000);
%      md=triangle(md,'DomainOutline.exp','Rifts.exp',1500);

%Figure out a characteristic area. Resolution is a node oriented concept (ex a 1000m  resolution node would 
%be made of 1000*1000 area squares). 
if (nargin==3),
	resolution=varargin{1};
	riftname='';
end
if (nargin==4),
	riftname=varargin{1};
	resolution=varargin{2};
end

%Check that mesh was not already run, and warn user: 
if md.mesh.numberofelements~=0,
	choice=input('This model already has a mesh. Are you sure you want to go ahead? (y/n)','s');
	if ~strcmp(choice,'y')
		disp('no meshing done ... exiting');
		return
	end
end

area=resolution^2;

%Check that file exist (this is a very very common mistake)
if ~exist(domainname)
	error(['file "' domainname '" not found']);
end

%Mesh using Triangle
[elements,x,y,segments,segmentmarkers]=Triangle_matlab(domainname,riftname,area);

%check that all the created nodes belong to at least one element
removeorphans=1;
if removeorphans,
	uniqueelements=sort(unique(elements(:)));
	orphans=find(~ismember([1:length(x)],uniqueelements));
	for i=1:length(orphans),
		disp('WARNING: removing orphans');
		%get rid of the orphan node i
		%update x and y
		x=[x(1:orphans(i)-(i-1)-1); x(orphans(i)-(i-1)+1:end)];
		y=[y(1:orphans(i)-(i-1)-1); y(orphans(i)-(i-1)+1:end)];
		%update elements
		pos=find(elements>orphans(i)-(i-1));
		elements(pos)=elements(pos)-1;
		%update segments
		pos1=find(segments(:,1)>orphans(i)-(i-1));
		pos2=find(segments(:,2)>orphans(i)-(i-1));
		segments(pos1,1)=segments(pos1,1)-1;
		segments(pos2,2)=segments(pos2,2)-1;
	end
end

%plug into md
md.mesh=mesh2d();
md.mesh.x=x;
md.mesh.y=y;
md.mesh.elements=elements;
md.mesh.segments=segments;
md.mesh.segmentmarkers=segmentmarkers;

%Fill in rest of fields:
md.mesh.numberofelements=size(md.mesh.elements,1);
md.mesh.numberofvertices=length(md.mesh.x);
md.mesh.vertexonboundary=zeros(md.mesh.numberofvertices,1);
md.mesh.vertexonboundary(md.mesh.segments(:,1:2))=1;

%Now, build the connectivity tables for this mesh.
md.mesh.vertexconnectivity=NodeConnectivity(md.mesh.elements,md.mesh.numberofvertices);
md.mesh.elementconnectivity=ElementConnectivity(md.mesh.elements,md.mesh.vertexconnectivity);
