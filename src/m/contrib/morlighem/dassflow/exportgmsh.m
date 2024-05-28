function exportgmsh(mesh,ocean_levelset,filename),
%EXPORTGMSH - export mesh to gmsh format
%
%   http://www.geuz.org/gmsh/doc/texinfo/#MSH-ASCII-file-format
%
%   Usage:
%      exportgmsh(mesh,ocean_levelset,filename)
%
%   Example:
%      exportgmsh(md.mesh,md.mask.ocean_levelset,'temp.msh')

fid=fopen(filename,'w');

%Header
fprintf(fid,'$MeshFormat\n');
fprintf(fid,'2.2 0 8\n');
fprintf(fid,'$EndMeshFormat\n');

%Vertices
nbv = mesh.numberofvertices;
fprintf(fid,'$Nodes\n');
fprintf(fid,'%i\n',nbv);
fprintf(fid,'%i %8g %8g %8g\n',[[1:nbv]' mesh.x mesh.y zeros(nbv,1)]');
fprintf(fid,'$EndNodes\n');

%Boundary Elements first
nbe     = mesh.numberofelements;
nbs     = size(mesh.segments,1);
segment = 1;
tria    = 2;

%Create flags
grounded = sum(ocean_levelset(mesh.segments(:,1:2))>0,2);
A = zeros(nbs,2);
pos = find(mesh.segmentmarkers==4);
A(pos,:)=repmat([1,4],[numel(pos) 1]);
pos = find(mesh.segmentmarkers==1 &  grounded);
A(pos,:)=repmat([2,1],[numel(pos) 1]);
pos = find(mesh.segmentmarkers==1 & ~grounded);
A(pos,:)=repmat([3,5],[numel(pos) 1]);
pos = find(mesh.segmentmarkers==2);
A(pos,:)=repmat([4,2],[numel(pos) 1]);
pos = find(mesh.segmentmarkers==3);
A(pos,:)=repmat([5,3],[numel(pos) 1]);

fprintf(fid,'$Elements\n');
fprintf(fid,'%i\n',nbe+nbs);
fprintf(fid,'%i %i %i %i %i %i %i\n',[[1    :nbs    ]' segment*ones(nbs,1) 2*ones(nbs,1) A mesh.segments(:,1:2)]');
fprintf(fid,'%i %i %i %i %i %i %i %i\n',[[nbs+1:nbs+nbe]' tria*ones(nbe,1) 2*ones(nbe,1) 7*ones(nbe,1) 6*ones(nbe,1) mesh.elements]');
fprintf(fid,'$EndElements\n');

fclose(fid);
