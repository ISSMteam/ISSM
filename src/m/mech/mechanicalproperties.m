function md=mechanicalproperties(md,vx,vy,varargin)
%MECHANICALPROPERTIES - compute stress and strain rate for a given velocity
%
%   this routine computes the components of the (deviatoric) stress tensor,
%   the strain rate tensor, and their respective principal directions.
%   The results are in the model md: md.results
%
%   Usage:
%      md=mechanicalproperties(md,vx,vy)
%
%   Example:
%      md=mechanicalproperties(md,md.initialization.vx,md.initialization.vy);
%      md=mechanicalproperties(md,md.inversion.vx_obs,md.inversion.vy_obs);

%some checks
if length(vx)~=md.mesh.numberofvertices | length(vy)~=md.mesh.numberofvertices,
	%error(['the input velocity should be of size ' num2str(md.mesh.numberofvertices) '!'])
end
if dimension(md.mesh)~=2
	error('only 2d model supported yet');
end
if any(md.flowequation.element_equation~=2),
	disp('Warning: the model has some non SSA elements. These will be treated like SSA''s elements');
end

%get damage, if passed
options = pairoptions(varargin{:});
damage = getfieldvalue(options,'damage',zeros(md.mesh.numberofvertices,1));

%initialization
numberofelements=md.mesh.numberofelements;
index=md.mesh.elements;
summation=[1;1;1];
directionsstress=zeros(numberofelements,4);
directionsstrain=zeros(numberofelements,4);
valuesstress=zeros(numberofelements,2);
valuesstrain=zeros(numberofelements,2);

%compute nodal functions coefficients N(x,y)=alpha x + beta y +gamma
[alpha beta]=GetNodalFunctionsCoeff(index,md.mesh.x,md.mesh.y);

%compute shear
vxlist=vx(index)/md.constants.yts;
vylist=vy(index)/md.constants.yts;
ux=(vxlist.*alpha)*summation;
uy=(vxlist.*beta)*summation;
vx=(vylist.*alpha)*summation;
vy=(vylist.*beta)*summation;						
uyvx=(vx+uy)./2;
clear vxlist vylist

%compute viscosity
nu=zeros(numberofelements,1);
B_bar=md.materials.rheology_B(index)*summation/3;
power=(md.materials.rheology_n-1)./(2*md.materials.rheology_n);
second_inv=(ux.^2+vy.^2+((uy+vx).^2)/4+ux.*vy);

%some corrections
location=find(second_inv==0 & power~=0);
nu(location)=10^18; 	%arbitrary maximum viscosity to apply where there is no effective shear

if isa(md.materials,'matice')
	location=find(second_inv~=0);
	nu(location)=B_bar(location)./(second_inv(location).^power(location));
	location=find(second_inv==0 & power==0);
	nu(location)=B_bar(location);
elseif isa(md.materials,'matdamageice')
	Zinv=1-damage(index)*summation/3;
	location=find(second_inv~=0);
	nu(location)=Zinv(location).*B_bar(location)./(second_inv(location).^power(location));
	location=find(second_inv==0 & power==0);
	nu(location)=Zinv(location).*B_bar(location);
	clear Zinv
else
	error(['class of md.materials (' class(md.materials) ') not recognized or not supported']);
end
clear B_bar location second_inv power

%compute stress
tau_xx=nu.*ux;
tau_yy=nu.*vy;
tau_xy=nu.*uyvx;

%compute principal properties of stress
for i=1:numberofelements,

	%compute stress and strainrate matrices
	stress=[tau_xx(i) tau_xy(i)
	tau_xy(i)  tau_yy(i)];
	strain=[ux(i) uyvx(i)
	uyvx(i)  vy(i)];

	%eigen values and vectors
	[directions,value]=eig(stress);
	%sort by algebraic value of eigenvalue (not absolute value) in descending order
	[val,idx]=sort(diag(value),'descend');
	%re-order eigenvalues and associated vectors 
	value=value(idx,idx);
	directions=directions(:,idx);
	valuesstress(i,:)=[value(1,1) value(2,2)];
	directionsstress(i,:)=directions(:)';
	[directions,value]=eig(strain);
	%same for strainrate
	[val,idx]=sort(diag(value),'descend');
	value=value(idx,idx);
	directions=directions(:,idx);
	valuesstrain(i,:)=[value(1,1) value(2,2)];
	directionsstrain(i,:)=directions(:)';
end

%plug onto the model
%NB: Matlab sorts the eigen value in increasing order, we want the reverse
strainrate=struct('xx',[],'yy',[],'xy',[],'principalvalue1',[],'principalaxis1',[],'principalvalue2',[],'principalaxis2',[],'effectivevalue',[]);
strainrate.xx=ux*md.constants.yts; %strain rate in 1/a instead of 1/s
strainrate.yy=vy*md.constants.yts; 
strainrate.xy=uyvx*md.constants.yts; 
strainrate.principalvalue1=valuesstrain(:,1)*md.constants.yts; 
strainrate.principalaxis1=directionsstrain(:,1:2);
strainrate.principalvalue2=valuesstrain(:,2)*md.constants.yts; 
strainrate.principalaxis2=directionsstrain(:,3:4);
strainrate.effectivevalue=1/sqrt(2)*sqrt(strainrate.xx.^2+strainrate.yy.^2+2*strainrate.xy.^2);
md.results.strainrate=strainrate;

%exact same stress as above
deviatoricstress=struct('xx',[],'yy',[],'xy',[],'principalvalue1',[],'principalaxis1',[],'principalvalue2',[],'principalaxis2',[],'effectivevalue',[]);
deviatoricstress.xx=tau_xx;
deviatoricstress.yy=tau_yy;
deviatoricstress.xy=tau_xy;
deviatoricstress.principalvalue1=valuesstress(:,1);
deviatoricstress.principalaxis1=directionsstress(:,1:2);
deviatoricstress.principalvalue2=valuesstress(:,2);
deviatoricstress.principalaxis2=directionsstress(:,3:4);
deviatoricstress.effectivevalue=1/sqrt(2)*sqrt(deviatoricstress.xx.^2+deviatoricstress.yy.^2+2*deviatoricstress.xy.^2);
md.results.deviatoricstress=deviatoricstress;

viscosity=struct('nu',[]);
viscosity.nu=nu;
md.results.viscosity=viscosity;
