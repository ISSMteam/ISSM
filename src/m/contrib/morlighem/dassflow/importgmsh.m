function md = importgmsh(filename,dim)

%some checks
if ~exist(filename),
	error(['expread error message: file ' filename ' not found!']);
end

%open file
fid=fopen(filename,'r');

%Get Mesh format
A=fscanf(fid,'%s',1);
if ~strcmp(A,'$MeshFormat'), 
	error(['Expecting $MeshFormat (' A ')']);
end
A=fscanf(fid,'%f %i %i',[1 3]);
A=fscanf(fid,'%s',1);
if ~strcmp(A,'$EndMeshFormat'), 
	error(['Expecting $EndMeshFormat (' A ')']);
end

%Nodes
A=fscanf(fid,'%s',1);
if ~strcmp(A,'$Nodes'), 
	error(['Expecting $Nodes (' A ')']);
end
nbv=fscanf(fid,'%i',1);
disp(['Number of nodes: ' num2str(nbv) ]);
A=fscanf(fid,'%i %f %f %f',[4 nbv]);
x = A(2,:)';
y = A(3,:)';
z = A(4,:)';

A=fscanf(fid,'%s',1);
if ~strcmp(A,'$EndNodes'), 
	error(['Expecting $EndNodes (' A ')']);
end

%Elements
A=fscanf(fid,'%s',1);
if ~strcmp(A,'$Elements'), 
	error(['Expecting $Elements (' A ')']);
end
nbt=fscanf(fid,'%i',1);
disp(['Number of elements: ' num2str(nbt) ]);
counter = 0;
if (dim==2),
	index   = zeros(0,3);
	segments       = zeros(0,2);
	segmentmarkers = zeros(0,1);
elseif (dim==3),
	index   = zeros(0,4);
	segments       = zeros(0,3);
	segmentmarkers = zeros(0,1);
else
	error('not supported');
end

while(counter<nbt);
	id = fscanf(fid,'%i',1);
	ty = fscanf(fid,'%i',1);
	nbf = fscanf(fid,'%i',1);
	flags = fscanf(fid,'%i',nbf);

	switch(ty)
		case 1, %segments
			A=fscanf(fid,'%i %i',2);
			if (dim==2),  %Actual element
				segments(end+1,:)=A;
				if    (flags(1)==5 & flags(2)==3), segmentmarkers(end+1)=3; 
				elseif(flags(1)==1 & flags(2)==4), segmentmarkers(end+1)=4;
				elseif(flags(1)==2 & flags(2)==1), segmentmarkers(end+1)=1;
				elseif(flags(1)==4 & flags(2)==2), segmentmarkers(end+1)=2;
				else error(['flags ' num2str(flags') ' not supported']);
				end
			else
				error('not supported');
			end
		case 2, %tria
			A=fscanf(fid,'%i %i %i',3);
			if (dim==2), %Actual element
				index(end+1,:)=A;
			else         %Boundary element
				segments(end+1,:)=A;
				if    (flags(1)==1), segmentmarkers(end+1)=1; 
				elseif(flags(1)==2), segmentmarkers(end+1)=2;
				elseif(flags(1)==3), segmentmarkers(end+1)=3;
				elseif(flags(1)==4), segmentmarkers(end+1)=4;
				else error(['flags ' num2str(flags') ' not supported']);
				end
			end
		case 4, %tetra
			A=fscanf(fid,'%i %i %i %i',4);
			if (dim==3), %Actual element
				index(end+1,:)=A;
			else
				error('not supported');
			end
		case 15, %point
			A=fscanf(fid,'%i',1);
			continue;
		otherwise,
			error(['Type ' num2str(ty) ' not supported']);
	end
	counter = counter + 1;
end

%recreate segments
if dim==2,
	nbs = size(segments,1);
	segments = [segments zeros(nbs,1)];
	for i=1:nbs,
		E = find(sum(ismember(index,segments(i,:)),2)>1);
		segments(i,3)=E;
	end
else
	nbs = size(segments,1);
	segments = [segments zeros(nbs,1)];
	for i=1:nbs,
		E = find(sum(ismember(index,segments(i,:)),2)>2);
		segments(i,4)=E;
	end
end

%close file
fclose(fid);

%Create model
if dim==2, %2d triangles
	md=meshconvert(model,index,x,y);
	md.mesh=mesh2dvertical(md.mesh);
	md.mesh.segmentmarkers=segmentmarkers;
	md.mesh.segments=segments;
	md.mesh.vertexonbase=zeros(md.mesh.numberofvertices,1);
	md.mesh.vertexonbase(find(vertexflags(md.mesh,1)))=1;
	md.mesh.vertexonsurface=zeros(md.mesh.numberofvertices,1);
	md.mesh.vertexonsurface(find(vertexflags(md.mesh,3)))=1;
else
	md=model();
	md.mesh=mesh3dtetras();
	md.mesh.x = x;
	md.mesh.y = y;
	md.mesh.z = z;
	md.mesh.elements = index;
	md.mesh.numberofelements=size(md.mesh.elements,1);
	md.mesh.numberofvertices=length(md.mesh.x);

	%base 2, surface 1, inflow 3, outflow 4
	md.mesh.vertexonbase=zeros(md.mesh.numberofvertices,1);
	md.mesh.vertexonbase(segments(find(segmentmarkers==2),1:3))=1;
	md.mesh.vertexonsurface=zeros(md.mesh.numberofvertices,1);
	md.mesh.vertexonsurface(segments(find(segmentmarkers==1),1:3))=1;
	md.mesh.vertexonboundary=zeros(md.mesh.numberofvertices,1);
	md.mesh.vertexonboundary(segments(:,1:3))=1;
end
