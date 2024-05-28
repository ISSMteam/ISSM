function dhdt=contourmassbalance(md,file)
%CONTOURMASSBALANCE - compute the mass balance on a contour
%
%   Usage:
%      dhdt=contourmassbalance(md,file)

%some checks
if nargin~=2,
	help contourmassbalance
	error('contourmassbalance error message: bad usage');
end
if ((length(md.initialization.vx)~=md.mesh.numberofvertices)|(length(md.initialization.vy)~=md.mesh.numberofvertices))
	error(['thicknessevolution error message: vx and vy should have a length of ' num2str(md.mesh.numberofvertices)])
end
if ~exist(file),
	error(['thicknessevolution error message: file ' file ' not found']);
end

%Get segments enveloping contour
segments=contourenvelope(md.mesh,file);
%md.stressbalance.icefront=segments; plotmodel(md,'data','pressureload','expdisp',file);

%get flag list of elements and nodes inside the contour
nodein=ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,file,'node',1);
elemin=(sum(nodein(md.mesh.elements),2)==size(md.mesh.elements,2));

%conputing Mass flux
x=md.mesh.x;
y=md.mesh.y;
vx=mean(md.initialization.vx(segments(:,1:end-1)),2);
vy=mean(md.initialization.vy(segments(:,1:end-1)),2);
H=mean(md.geometry.thickness(segments(:,1:end-1)),2);
nx=cos(atan2((x(segments(:,1))-x(segments(:,2))) , (y(segments(:,2))-y(segments(:,1)))));
ny=sin(atan2((x(segments(:,1))-x(segments(:,2))) , (y(segments(:,2))-y(segments(:,1)))));
L=sqrt((x(segments(:,1))-x(segments(:,2))).^2+(y(segments(:,2))-y(segments(:,1))).^2);
flux = - md.materials.rho_ice*sum(L.*H.*(vx.*nx+vy.*ny)); %outflux is negative!
disp(['mass outflux on ' file ' = ' num2str(-flux/10^9) ' Gt/yr']);
areas=GetAreas(md.mesh.elements,md.mesh.x,md.mesh.y);
dhdt=flux/(sum(areas(find(elemin)))*md.materials.rho_ice);
disp(['dhdt on ' file ' (Flux  method) = ' num2str(dhdt) ' m/yr']);

dhdt=thicknessevolution(md);
in=find(elemin);
dhdt=sum(dhdt(in).*areas(in))/sum(areas(in));
disp(['dhdt on ' file ' (divHV method) = ' num2str(dhdt) ' m/yr']);
