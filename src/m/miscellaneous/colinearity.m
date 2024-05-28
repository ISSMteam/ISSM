function theta = colinearity(md)

%load some variables (it is much faster if the variab;es are loaded from md once for all) 
if ~strcmpi(meshtype(md.mesh),'3D'),
	numberofelements=md.mesh.numberofelements;
	numberofnodes=md.mesh.numberofvertices;
	index=md.mesh.elements;
	x=md.mesh.x; y=md.mesh.y;
else
	numberofelements=md.mesh.numberofelements2d;
	numberofnodes=md.mesh.numberofvertices2d;
	index=md.mesh.elements2d;
	x=md.mesh.x2d; y=md.mesh.y2d;
end

%compute nodal functions coefficients N(x,y)=alpha x + beta y + gamma
[alpha beta]=GetNodalFunctionsCoeff(index,x,y);

s = averaging(md,md.geometry.surface,2);

summation=[1;1;1];
dsdx=(s(index).*alpha)*summation;
dsdy=(s(index).*beta)*summation;
dsdx = -averaging(md,dsdx,0);
dsdy = -averaging(md,dsdy,0);

vx = md.inversion.vx_obs;
vy = md.inversion.vy_obs;
v  = md.inversion.vel_obs;
v2 = sqrt(dsdx.^2 + dsdy.^2);

theta = acos((vx.*dsdx + vy.*dsdy)./(v.*v2+eps));
