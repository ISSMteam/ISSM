function md=BamgCall(md,field,hmin,hmax,gradation,epsilon),
%BAMGCALL - call bam
%
%   build a metric using the Hessian of the given field
%   call Bamg and the output mesh is plugged onto the model
%   -hmin = minimum edge length (m)
%   -hmax = maximum edge length (m)
%   -gradation = maximum edge length gradation between 2 elements
%   -epsilon = average error on each element (m/yr)
%
%   Usage:
%      md=BamgCall(md,field,hmin,hmax,gradation,epsilon);
%
%   Example:
%      md=BamgCall(md,md.inversion.vel_obs,1500,10^8,1.3,0.9);

%2d geometric parameter (do not change)
scale=2/9; 

%Compute Hessian
t1=clock; fprintf('%s','      computing Hessian...');
hessian=ComputeHessian(md.mesh.elements,md.mesh.x,md.mesh.y,field,'node');
t2=clock;fprintf('%s\n',[' done (' num2str(etime(t2,t1)) ' seconds)']);

%Compute metric
t1=clock; fprintf('%s','      computing metric...');
metric=ComputeMetric(hessian,scale,epsilon,hmin,hmax,[]);
t2=clock;fprintf('%s\n',[' done (' num2str(etime(t2,t1)) ' seconds)']);

%write files
t1=clock; fprintf('%s','      writing initial mesh files...');
fid=fopen('carre0.met','w');
fprintf(fid,'%i %i\n',md.mesh.numberofvertices,3);
fprintf(fid,'%i %i %i\n',metric');
fclose(fid);

fid=fopen('carre0.mesh','w');

%initialiation
fprintf(fid,'%s %i\n','MeshVersionFormatted',0);

%dimension
fprintf(fid,'\n%s\n%i\n','Dimension',2);

%Vertices
fprintf(fid,'\n%s\n%i\n\n','Vertices',md.mesh.numberofvertices);
fprintf(fid,'%8g %8g %i\n',[md.mesh.x md.mesh.y ones(md.mesh.numberofvertices,1)]');

%Triangles
fprintf(fid,'\n\n%s\n%i\n\n','Triangles',md.mesh.numberofelements);
fprintf(fid,'%i %i %i %i\n',[md.mesh.elements ones(md.mesh.numberofelements,1)]');
numberofelements1=md.mesh.numberofelements;

%close
fclose(fid);
t2=clock;fprintf('%s\n',[' done (' num2str(etime(t2,t1)) ' seconds)']);

%call bamg
fprintf('%s\n','      call Bamg...');
system(['bamg -ratio ' num2str(gradation) ' -splitpbedge -nbv 1000000 -M carre0.met -b carre0.mesh -o carre1.mesh']);

%plug new mesh
t1=clock; fprintf('\n%s','      reading final mesh files...');
A=meshread('carre1.mesh');
md.mesh.x=A.x;
md.mesh.y=A.y;
md.z=zeros(A.nods,1);
md.mesh.elements=A.index;
md.mesh.numberofvertices=A.nods;
md.mesh.numberofelements=A.nels;
numberofelements2=md.mesh.numberofelements;
t2=clock;fprintf('%s\n\n',[' done (' num2str(etime(t2,t1)) ' seconds)']);

%display number of elements
fprintf('\n%s %i','      inital number of elements:',numberofelements1);
fprintf('\n%s %i\n\n','      new    number of elements:',numberofelements2);

%clean up:
system('rm carre0.mesh carre0.met carre1.mesh carre1.mesh.gmsh');
