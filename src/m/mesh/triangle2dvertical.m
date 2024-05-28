function md=triangle(md,domainname,resolution)
%TRIANGLE - create model mesh using the triangle package
%
%   This routine creates a model mesh using Triangle and a domain outline, to within a certain resolution
%   where md is a @model object, domainname is the name of an Argus domain outline file, 
%   and resolution is a characteristic length for the mesh (same unit as the domain outline
%   unit). Riftname is an optional argument (Argus domain outline) describing rifts.
%
%   Usage:
%      md=triangle(md,domainname,resolution)
%
%   Examples:
%      md=triangle(md,'DomainOutline.exp',1000);

%Check that mesh was not already run, and warn user: 
if md.mesh.numberofelements~=0,
	choice=input('This model already has a mesh. Are you sure you want to go ahead? (y/n)','s');
	if ~strcmp(choice,'y')
		disp('no meshing done ... exiting');
		return
	end
end

area=resolution^2;

%Mesh using Triangle
[elements,x,z,segments,segmentmarkers]=Triangle_matlab(domainname,'',area);

%check that all the created nodes belong to at least one element
orphan=find(~ismember([1:length(x)],sort(unique(elements(:)))));
for i=1:length(orphan),
	disp('WARNING: removing orphans');
	%get rid of the orphan node i
	%update x and y
	x=[x(1:orphan(i)-(i-1)-1); x(orphan(i)-(i-1)+1:end)];
	z=[z(1:orphan(i)-(i-1)-1); z(orphan(i)-(i-1)+1:end)];
	%update elements
	pos=find(elements>orphan(i)-(i-1));
	elements(pos)=elements(pos)-1;
	%update segments
	pos1=find(segments(:,1)>orphan(i)-(i-1));
	pos2=find(segments(:,2)>orphan(i)-(i-1));
	segments(pos1,1)=segments(pos1,1)-1;
	segments(pos2,2)=segments(pos2,2)-1;
end

%plug into md
md.mesh=mesh2dvertical();
md.mesh.x=x;
md.mesh.z=z;
md.mesh.elements=elements;
md.mesh.segments=segments;
md.mesh.segmentmarkers=segmentmarkers;

%Fill in rest of fields:
md.mesh.numberofelements=size(md.mesh.elements,1);
md.mesh.numberofvertices=length(md.mesh.x);
md.mesh.vertexonboundary=zeros(md.mesh.numberofvertices,1); md.mesh.vertexonboundary(md.mesh.segments(:,1:2))=1;

%Now, build the connectivity tables for this mesh.
md.mesh.vertexconnectivity=NodeConnectivity(md.mesh.elements,md.mesh.numberofvertices);
md.mesh.elementconnectivity=ElementConnectivity(md.mesh.elements,md.mesh.vertexconnectivity);
