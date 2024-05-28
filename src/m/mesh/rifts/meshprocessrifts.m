function md=meshprocessrifts(md,domainoutline)
%MESHPROCESSRIFTS - process mesh when rifts are present
%
%   split rifts inside mesh (rifts are defined by presence of
%   segments inside the domain outline)
%   if domain outline is provided, check for rifts that could touch it, and open them up.
%
%   Usage:
%      md=meshprocessrifts(md,domainoutline)
%
%   Ex: 
%      md=meshprocessrifts(md,'DomainOutline.exp');
%

%some checks on arguments: 
if nargout~=1,
	help meshprocessrifts
	error('meshprocessrifts usage error:');
end

if nargin~=2,
	help meshprocessrifts
	error('meshprocessrifts usage error:');
end

%Call MEX file
[md.mesh.elements,md.mesh.x,md.mesh.y,md.mesh.segments,md.mesh.segmentmarkers,md.rifts.riftstruct]=ProcessRifts(md.mesh.elements,md.mesh.x,md.mesh.y,md.mesh.segments,md.mesh.segmentmarkers);
if ~isstruct(md.rifts.riftstruct),
	error('ProcessRifts did not find any rift');
end

%Fill in rest of fields:
numrifts=length(md.rifts.riftstruct);
md.mesh.numberofelements=length(md.mesh.elements);
md.mesh.numberofvertices=length(md.mesh.x);
md.mesh.vertexonboundary=zeros(length(md.mesh.x),1); md.mesh.vertexonboundary(md.mesh.segments(:,1:2))=1;

%get coordinates of rift tips
for i=1:numrifts,
	md.rifts.riftstruct(i).tip1coordinates=[md.mesh.x(md.rifts.riftstruct(i).tips(1)) md.mesh.y(md.rifts.riftstruct(i).tips(1))];
	md.rifts.riftstruct(i).tip2coordinates=[md.mesh.x(md.rifts.riftstruct(i).tips(2)) md.mesh.y(md.rifts.riftstruct(i).tips(2))];
end

%In case we have rifts that open up the domain outline, we need to open them: 
flags=ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,domainoutline,'node',0);
found=0;
for i=1:numrifts,
	if flags(md.rifts.riftstruct(i).tips(1))==0,
		found=1;
		break;
	end
	if flags(md.rifts.riftstruct(i).tips(2))==0,
		found=1;
		break;
	end
end
if found,
	md=meshprocessoutsiderifts(md,domainoutline);
end

%get elements that are not correctly oriented in the correct direction:
aires=GetAreas(md.mesh.elements,md.mesh.x,md.mesh.y);
pos=find(aires<0);
md.mesh.elements(pos,:)=[md.mesh.elements(pos,2) md.mesh.elements(pos,1) md.mesh.elements(pos,3)];

%case of 3D surface mesh: 
if strcmpi(class(md.mesh),'mesh3dsurface'),
	md.mesh.z=md.mesh.x; md.mesh.z(:)=0;
end
