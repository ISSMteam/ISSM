function ExportGmsh(md,filename)
%EXPORTGMSH - export mesh to gmsh format
%
%   Usage:
%      ExportGmsh(md,filename)

t1=clock;fprintf('%s',['writing gmsh mesh file']);
fid=fopen(filename,'w');

%initialization
fprintf(fid,'$MeshFormat \n');
fprintf(fid,'2.2 0 8 \n');
fprintf(fid,'$EndMeshFormat \n');
fprintf(fid,'$Nodes \n');
fprintf(fid,'%i \n',md.mesh.numberofvertices);
np=0;
%printing point positions
for ii=1:md.mesh.numberofvertices
	np = np+1;
	fprintf(fid,'%g %14.7e %14.7e 0.0 \n',np,md.mesh.x(np),md.mesh.y(np));
end

fprintf(fid,'$EndNodes \n');
fprintf(fid,'$Elements \n');
fprintf(fid,'%i \n',md.mesh.numberofelements+size(md.mesh.segments,1));
np=0;

%printing elements caracteristics for boundaries

for ii=1:size(md.mesh.segments,1)
	np = np+1;
	if(md.mesh.x(md.mesh.segments(np,1))==max(md.mesh.x(:))&&md.mesh.x(md.mesh.segments(np,2))==max(md.mesh.x(:))),
		bc_id=1;
	elseif(md.mesh.y(md.mesh.segments(np,1))==max(md.mesh.y(:))&&md.mesh.y(md.mesh.segments(np,2))==max(md.mesh.y(:))),
		bc_id=2;
	elseif(md.mesh.x(md.mesh.segments(np,1))==min(md.mesh.x(:))&&md.mesh.x(md.mesh.segments(np,2))==min(md.mesh.x(:))),
		bc_id=3;
	elseif(md.mesh.y(md.mesh.segments(np,1))==min(md.mesh.y(:))&&md.mesh.y(md.mesh.segments(np,2))==min(md.mesh.y(:))),
		bc_id=4;
	else
		bc_id=5;
  end
  fprintf(fid,'%g 1 2 %g 1 %g %g \n',np,bc_id,md.mesh.segments(np,1),md.mesh.segments(np,2));
end
%and for the body
body_id=1;
for ii=1:md.mesh.numberofelements
  np = np+1;
  fprintf(fid,'%g 2 2 %g 3 %g %g %g \n',np,body_id,md.mesh.elements(ii,1),md.mesh.elements(ii,2),md.mesh.elements(ii,3));
end
fprintf(fid,'$EndElements \n');
%close
fclose(fid);
t2=clock;fprintf('%s\n',[' done (' num2str(etime(t2,t1)) ' seconds)']);
