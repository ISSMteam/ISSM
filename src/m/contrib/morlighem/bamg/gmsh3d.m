function md=gmsh3d(md,domainfile,h),

%Read domain
domain=expread(domainfile);
x   = domain.x(1:end-1);
y   = domain.y(1:end-1);
nbv = numel(x);

%write files
t1=clock; fprintf('%s','      writing initial mesh files...');
fid=fopen('model.geo','w');
fprintf(fid,['// Gmsh input file, created by ISSM on ' date '\n']);

%Lower points
for i=1:nbv,
	fprintf(fid,'Point(%i) = {%8g, %8g, %8g};\n',i,x(i),y(i),0.);
end
%Upper points
for i=1:nbv,
	fprintf(fid,'Point(%i) = {%8g, %8g, %8g};\n',nbv+i,x(i),y(i),1.);
end

%Lower lines
for i=1:nbv-1
	fprintf(fid,'Line(%i) = {%i, %i};\n',i,i,i+1);
end
fprintf(fid,'Line(%i) = {%i, %i};\n',nbv,nbv,1);
%Upper lines
for i=nbv+1:2*nbv-1
	fprintf(fid,'Line(%i) = {%i, %i};\n',i,i,i+1);
end
fprintf(fid,'Line(%i) = {%i, %i};\n',2*nbv,2*nbv,nbv+1);
%Side lines
for i=1:nbv
	fprintf(fid,'Line(%i) = {%i, %i};\n',2*nbv+i,i,nbv+i);
end

counter = 3*nbv;
ps = zeros(nbv+2,1);
%Lower surface Loop and surface
counter = counter+1;
fprintf(fid,['Line Loop(' num2str(counter) ') = {']);
for i=1:nbv-1
	fprintf(fid,'%i,',i);
end
fprintf(fid,'%i};\n',nbv);
fprintf(fid,['Plane Surface(' num2str(counter+1) ') = {' num2str(counter) '};\n']);
ps(1)=counter+1;
counter = counter+1;
%Upper surface Loop and surface
counter = counter+1;
fprintf(fid,['Line Loop(' num2str(counter) ') = {']);
for i=nbv+1:2*nbv-1
	fprintf(fid,'%i,',i);
end
fprintf(fid,'%i};\n',2*nbv);
ps(2)=counter+1;
fprintf(fid,['Plane Surface(' num2str(counter+1) ') = {' num2str(counter) '};\n']);
counter = counter+2;
%Sides surfaces
for i=1:nbv-1,
	fprintf(fid,['Line Loop(' num2str(counter) ') = {' num2str(i) ',' num2str(2*nbv+i+1) ',-' num2str(nbv+i) ',-' num2str(2*nbv+i) '};\n']);
	fprintf(fid,['Plane Surface(' num2str(counter+1) ') = {' num2str(counter) '};\n']);
	ps(2+i)=counter+1;
	counter=counter+2;
end
fprintf(fid,['Line Loop(' num2str(counter) ') = {' num2str(nbv) ',' num2str(2*nbv+1) ',-' num2str(2*nbv) ',-' num2str(3*nbv) '};\n']);
fprintf(fid,['Plane Surface(' num2str(counter+1) ') = {' num2str(counter) '};\n']);
ps(2+nbv)=counter+1;
counter=counter+2;

%Physical surfaces
counter = counter+1;
fprintf(fid,['Surface Loop(' num2str(counter) ') = {']);
for i=1:numel(ps)-1
	fprintf(fid,'%i,',ps(i));
end
fprintf(fid,'%i};\n',ps(end));
fprintf(fid,['Physical Surface(1) = {' num2str(ps(2)) '};\n']);
fprintf(fid,['Physical Surface(2) = {' num2str(ps(1)) '};\n']);
fprintf(fid,['Physical Surface(3) = {']);
for i=3:numel(ps)-1
	fprintf(fid,'%i,',ps(i));
end
fprintf(fid,'%i};\n',ps(end));

%Volume
fprintf(fid,['Volume(1) = {' num2str(counter) '};\n']);
fprintf(fid,['Physical Volume(2) = {1};\n']);

%resolution
fprintf(fid,'Mesh.CharacteristicLengthMax = %g;',h);
fclose(fid);

t2=clock;fprintf('%s\n',[' done (' num2str(etime(t2,t1)) ' seconds)']);

%call gmsh
fprintf('%s\n','      call gmsh...');
system([issmdir() '/externalpackages/gmsh/install/gmsh -3 -v 0 model.geo']);

%plug new mesh
t1=clock; fprintf('\n%s','      reading final mesh files...');
md=importgmsh('model.msh',3);
