function md=gmsh(md,domainfile,h),

%Read domain
domain=expread(domainfile);
x   = domain.x(1:end-1);
y   = domain.y(1:end-1);
nbv = numel(x);

%write files
t1=clock; fprintf('%s','      writing initial mesh files...');
fid=fopen('model.geo','w');
fprintf(fid,['// Gmsh input file, created by ISSM on ' date '\n']);
for i=1:nbv,
	fprintf(fid,'Point(%i) = {%8g, %8g, %8g};\n',i,x(i),y(i),0.);
end
for i=1:nbv-1
	fprintf(fid,'Line(%i) = {%i, %i};\n',i,i,i+1);
end
fprintf(fid,'Line(%i) = {%i, %i};\n',nbv,nbv,1);
fprintf(fid,'Line Loop(5) = {');
for i=1:nbv-1
	fprintf(fid,'%i,',i);
end
fprintf(fid,'%i};\n',nbv);
fprintf(fid,'Plane Surface(6) = {5};\n');

%Physical lines and surfaces
fprintf(fid,'Physical Line(2) = {1};\n');
fprintf(fid,'Physical Line(4) = {2};\n');
fprintf(fid,'Physical Line(5) = {3};\n');
fprintf(fid,'Physical Line(1) = {4};\n');
fprintf(fid,'Physical Surface(7) = {6};\n');

%resolution
fprintf(fid,'Mesh.CharacteristicLengthMax = %g;',h);

%fprintf(fid,'Plane Surface(7) = {6, 2};\n');
fclose(fid);

t2=clock;fprintf('%s\n',[' done (' num2str(etime(t2,t1)) ' seconds)']);

%call gmsh
fprintf('%s\n','      call gmsh...');
system([issmdir() '/externalpackages/gmsh/install/gmsh -2 model.geo']);

%plug new mesh
t1=clock; fprintf('\n%s','      reading final mesh files...');
md=importgmsh('model.msh',2);
