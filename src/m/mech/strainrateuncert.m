function md=strainratuncert(md,vx,vy,dvx,dvy)
%STRAINRATEUNCERT - compute uncertainty in strain rate components
%
%   this routine computes the uncertainties in the strain rate tensor
%	 components given the uncertainty in surface velocity data.
%   The results are stored in md.results
%
%	 'dvx' and 'dvy' are velocity errors in x and y components in m/yr.  
%   These can either be scalars or arrays of length md.mesh.numberofvertices
%
%   Usage:
%      md=strainrateuncert(md,vx,vy,dv)
%
%   Example:
%      md=mechanicalproperties(md,md.initialization.vx,md.initialization.vy,5);
%      md=mechanicalproperties(md,md.inversion.vx_obs,md.inversion.vy_obs,dv);

%some checks
if length(vx)~=md.mesh.numberofvertices | length(vy)~=md.mesh.numberofvertices,
	error(['the input velocity should be of size ' num2str(md.mesh.numberofvertices) '!'])
end
if length(dvx)==1,
	dvx=dvx*ones(md.mesh.numberofelements,1);
end
if length(dvx)~=md.mesh.numberofelements,
	error(['the velocity error dvx should be of size ' num2str(md.mesh.numberofelements) ' or 1!'])
end
if length(dvy)==1,
	dvy=dvy*ones(md.mesh.numberofelements,1);
end
if length(dvy)~=md.mesh.numberofelements,
	error(['the velocity error dvy should be of size ' num2str(md.mesh.numberofelements) ' or 1!'])
end
if dimension(md.mesh)~=2,
	error('only 2d model supported yet');
end
if any(md.flowequation.element_equation~=2),
	disp('Warning: the model has some non SSA elements. These will be treated like SSA''s elements');
end

%initialization
index=md.mesh.elements;
summation=[1;1;1];
dvxlist=dvx(index);
dvylist=dvy(index);

%compute nodal functions coefficients N(x,y)=alpha x + beta y +gamma
[alpha beta]=GetNodalFunctionsCoeff(index,md.mesh.x,md.mesh.y);

strainrateuncert=struct('xx',[],'yy',[],'xy',[],'principalvalue1',[],'principalvalue2',[],'effectivevalue',[]);

strainrateuncert.xx=sqrt((dvxlist.*alpha).^2*summation);
strainrateuncert.yy=sqrt((dvylist.*beta).^2*summation);
strainrateuncert.xy=0.5*sqrt((dvxlist.*beta).^2*summation+(dvylist.*alpha).^2*summation);

exx=md.results.strainrate.xx; 
eyy=md.results.strainrate.yy;
exy=md.results.strainrate.xy;
p1a=strainrateuncert.xx.*(0.5+0.25*(((exx-eyy)/2).^2+exy.^2).^(-1./2).*(exx-eyy));
p2a=strainrateuncert.yy.*(0.5-0.25*(((exx-eyy)/2).^2+exy.^2).^(-1./2).*(exx-eyy));
p3a=strainrateuncert.xy.*((((exx-eyy)/2).^2+exy.^2).^(-1./2).*exy);
p1b=strainrateuncert.xx.*(0.5-0.25*(((exx-eyy)/2).^2+exy.^2).^(-1./2).*(exx-eyy));
p2b=strainrateuncert.yy.*(0.5+0.25*(((exx-eyy)/2).^2+exy.^2).^(-1./2).*(exx-eyy));
p3b=strainrateuncert.xy.*(-(((exx-eyy)/2).^2+exy.^2).^(-1./2).*exy);
strainrateuncert.principalvalue1=sqrt(p1a.^2+p2a.^2+p3a.^2);
strainrateuncert.principalvalue2=sqrt(p1b.^2+p2b.^2+p3b.^2);

effa=strainrateuncert.xx/sqrt(2).*(exx.^2+eyy.^2+2*exy.^2).^(-1./2).*exx;
effb=strainrateuncert.yy/sqrt(2).*(exx.^2+eyy.^2+2*exy.^2).^(-1./2).*eyy;
effc=2*strainrateuncert.xy/sqrt(2).*(exx.^2+eyy.^2+2*exy.^2).^(-1./2).*exy;
strainrateuncert.effectivevalue=sqrt(effa.^2+effb.^2+effc.^2);

md.results.strainrateuncert=strainrateuncert;
