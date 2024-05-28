function md=YamsCall(md,field,hmin,hmax,gradation,epsilon),
%YAMSCALL - call yams
%
%   build a metric using the Hessian of the given field
%   call Yams and the output mesh is plugged onto the model
%   -hmin = minimum edge length (m)
%   -hmax = maximum edge length (m)
%   -gradation = maximum edge length gradation between 2 elements
%   -epsilon = average error on each element (m/yr)
%
%   Usage:
%      md=YamsCall(md,field,hmin,hmax,gradation,epsilon);
%
%   Example:
%      md=YamsCall(md,md.inversion.vel_obs,1500,10^8,1.3,0.9);

%2d geometric parameter (do not change)
scale=2./9.;

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
save -ascii carre0.met  metric

fid=fopen('carre0.mesh','w');

%initialiation
fprintf(fid,'\n%s\n%i\n','MeshVersionFormatted',1);

%dimension
fprintf(fid,'\n%s\n%i\n','Dimension',2);

%Vertices
fprintf(fid,'\n%s\n%i\n\n','Vertices',md.mesh.numberofvertices);
fprintf(fid,'%8g %8g %i\n',[md.mesh.x md.mesh.y zeros(md.mesh.numberofvertices,1)]');

%Triangles
fprintf(fid,'\n\n%s\n%i\n\n','Triangles',md.mesh.numberofelements);
fprintf(fid,'%i %i %i %i\n',[md.mesh.elements zeros(md.mesh.numberofelements,1)]');
numberofelements1=md.mesh.numberofelements;

%Deal with rifts
if ~isnan(md.rifts.riftstruct),

	%we have the list of triangles that make up the rift. keep those triangles around during refinement.
	triangles=[];
	for i=1:size(md.rifts.riftstruct,1),
		triangles=[triangles md.rifts(i).riftstruct.segments(:,3)'];
	end

	fprintf(fid,'\n\n%s\n%i\n\n','RequiredTriangles',length(triangles));
	fprintf(fid,'%i\n',triangles);
end

%close
fclose(fid);
t2=clock;fprintf('%s\n',[' done (' num2str(etime(t2,t1)) ' seconds)']);

%call yams
fprintf('%s\n','      call Yams...');
if ispc()
	%windows
	system(['yams2-win -O 1 -v -0 -ecp -hgrad ' num2str(gradation)  ' carre0 carre1']);
elseif ismac()
	%Macosx
	system(['yams2-osx -O 1 -v -0 -ecp -hgrad ' num2str(gradation)  ' carre0 carre1']);
else
	%Linux
	system(['yams2-linux -O 1 -v -0 -ecp -hgrad ' num2str(gradation)  ' carre0 carre1']);
end

%plug new mesh
t1=clock; fprintf('\n%s','      reading final mesh files...');
Tria=load('carre1.tria');
Coor=load('carre1.coor');
md.mesh.x=Coor(:,1);
md.mesh.y=Coor(:,2);
md.mesh.elements=Tria;
md.mesh.numberofvertices=size(Coor,1);
md.mesh.numberofelements=size(Tria,1);
numberofelements2=md.mesh.numberofelements;
t2=clock;fprintf('%s\n\n',[' done (' num2str(etime(t2,t1)) ' seconds)']);

%display number of elements
fprintf('\n%s %i','      inital number of elements:',numberofelements1);
fprintf('\n%s %i\n\n','      new    number of elements:',numberofelements2);

%clean up:
system('rm carre0.mesh carre0.met carre1.tria carre1.coor carre1.meshb');
