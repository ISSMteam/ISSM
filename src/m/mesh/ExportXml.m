function ExportXml(md,filename)
%EXPORTGMSH - export mesh to xml format (For FEniCS)
%
%   Usage:
%      ExportXml(md,filename)

t1=clock;fprintf('%s',['writing xml mesh file']);
fid=fopen(filename,'w');

%initialization
fprintf(fid,'<?xml version="1.0" encoding="UTF-8"?>\n');
fprintf(fid,'\n');
fprintf(fid,'<dolfin xmlns:dolfin="http://www.fenicsproject.org">\n');
if isa(md.mesh,'mesh2d')
	fprintf(fid,'  <mesh celltype="triangle" dim="2">\n');
elseif isa(md.mesh,'mesh3dtetras')
	fprintf(fid,'  <mesh celltype="tetrahedron" dim="3">\n');
else
	error('only triangles and tets are supported');
end

%printing point positions
fprintf(fid,'    <vertices size="%i">\n',md.mesh.numberofvertices);
for i=1:md.mesh.numberofvertices
	fprintf(fid,'      <vertex index="%i" x="%17.15e" y="%17.15e" z="%17.15e"/>\n',i-1,md.mesh.x(i),md.mesh.y(i),md.mesh.z(i));
end
fprintf(fid,'    </vertices>\n');
fprintf(fid,'    <cells size="%i">\n',md.mesh.numberofelements);
if isa(md.mesh,'mesh2d')
	for i=1:md.mesh.numberofelements
		fprintf(fid,'      <triangle index="%i" v0="%i" v1="%i" v2="%i"/>\n',i-1,md.mesh.elements(i,1)-1,md.mesh.elements(i,2)-1,md.mesh.elements(i,3)-1);
	end
elseif isa(md.mesh,'mesh3dtetras')
	for i=1:md.mesh.numberofelements
		fprintf(fid,'      <tetrahedron index="%i" v0="%i" v1="%i" v2="%i" v3="%i"/>\n',i-1,md.mesh.elements(i,1)-1,md.mesh.elements(i,2)-1,md.mesh.elements(i,3)-1,md.mesh.elements(i,4)-1);
	end
else
	error('only triangles and tets are supported');
end
fprintf(fid,'    </cells>\n');
fprintf(fid,'  </mesh>\n');
fprintf(fid,'</dolfin>\n');

%close
fclose(fid);
t2=clock;fprintf('%s\n',[' done (' num2str(etime(t2,t1)) ' seconds)']);
