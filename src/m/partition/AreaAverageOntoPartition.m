function partvector=AreaAverageOntoPartition(md,vector,partition,layer)
%AREAAVERAGEONTOPARTITION 
%   compute partition values for a certain vector expressed on the vertices of the mesh.
%   Use area weighted average.
%
%   Usage:
%      average=AreaAverageOntoPartition(md,vector)
%      average=AreaAverageOntoPartition(md,vector,layer) %if in 3D, chose which layer is partitioned

%some checks
if dimension(md.mesh)==3,
	if nargin~=3,
		error('layer should be provided onto which Area Averaging occurs');
	end
	%save 3D model
	md3d=md;

	md.mesh.elements=md.mesh.elements2d;
	md.mesh.x=md.mesh.x2d;
	md.mesh.y=md.mesh.y2d;
	md.mesh.numberofvertices=md.mesh.numberofvertices2d;
	md.mesh.numberofelements=md.mesh.numberofelements2d;
	md.qmu.vertex_weight=[];
	md.mesh.vertexconnectivity=[];

	%run connectivity routine
	md=adjacency(md);

	%finally, project vector: 
	vector=project2d(md3d,vector,layer);
	partition=project2d(md3d,partition,layer);
end

%ok, first check that part is Matlab indexed
part=partition+1;

%some check: 
npart=qmupart2npart(partition);
if npart~=max(part),
	error('AreaAverageOntoPartition error message: ''npart'' should be equal to max(partition)');
end

%initialize output
partvector=zeros(max(part),1);

%start weight average
weightedvector=vector.*md.qmu.vertex_weight;
for i=1:max(part),
	pos=find(part==i);
	partvector(i)=sum(weightedvector(pos))/sum(md.qmu.vertex_weight(pos));
end

%in 3D, restore 3D model:
if dimension(md.mesh)==3,
	md=md3d;
end
