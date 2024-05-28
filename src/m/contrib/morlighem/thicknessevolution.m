function dhdt=thicknessevolution(md)
%THICKNESSEVOLUTION - compute the new thickness of a model after ∆t
%
%   This routine compute the new thickness of a model after a time step
%   according to the following formula:
%   dh/dt=-div(Hu)
%
%   Usage:
%      dhdt=thicknessevolution(md)

if (length(md.initialization.vx)~=md.mesh.numberofvertices)|(length(md.initialization.vy)~=md.mesh.numberofvertices)
	error(['thicknessevolution error message: vx and vy should have a length of ' num2str(md.mesh.numberofvertices)])
end

%load some variables 
H=md.geometry.thickness;
vx=md.initialization.vx;
vy=md.initialization.vy;
index=md.mesh.elements;

%compute nodal functions coefficients N(x,y)=alpha x + beta y + gamma
[alpha beta]=GetNodalFunctionsCoeff(md.mesh.elements,md.mesh.x,md.mesh.y); 

%compute dhdt=div(Hu)
summation=1/3*ones(3,1);
dhdt=(vx(index)*summation).*sum( H(index).*alpha,2) + (vy(index)*summation).*sum(H(index).*beta,2) ...
	+ ( H(index)*summation).*sum(vx(index).*alpha,2) + ( H(index)*summation).*sum(vy(index).*beta,2);
dhdt=-dhdt;
