function [alpha beta varargout]=GetNodalFunctionsCoeff(index,x,y)
%GETNODELFUNCTIONSCOEFF - compute nodal functions coefficients
%
%   Compute the coefficients alpha beta and optionaly gamma of
%   2d triangular elements. For each element, the nodal function
%   is defined as:
%   N(x,y)=sum(i=1:3) alpha_i * x + beta_i * y + gamma_i
%
%   Usage:
%      [alpha beta]=GetNodalFunctionsCoeff(index,x,y);
%      [alpha beta gamma]=GetNodalFunctionsCoeff(index,x,y);
%
%   Example:
%      [alpha beta gamma]=GetNodalFunctionsCoeff(md.mesh.elements,md.mesh.x,md.mesh.y);

%make columns out of x and y
x=x(:); y=y(:);

%get nels and nods
nels=size(index,1);
nods=length(x);

%some checks
if nargin~=3 | (nargout~=2 & nargout~=3),
	help GetNodalFunctionsCoeff
	error('GetNodalFunctionsCoeff error message: bad usage')
end
if length(y)~=nods,
	error('GetNodalFunctionsCoeff error message: x and y do not have the same length')
end
if max(index(:))>nods,
	error(['GetNodalFunctionsCoeff error message: index should not have values above ' num2str(nods) ])
end
if size(index,2)~=3,
	error('GetNodalFunctionsCoeff error message: only 2d meshes supported. index should have 3 columns.')
end

%initialize output
alpha=zeros(nels,3);
beta=zeros(nels,3);

%compute nodal functions coefficients N(x,y)=alpha x + beta y +gamma
x1=x(index(:,1)); x2=x(index(:,2)); x3=x(index(:,3));
y1=y(index(:,1)); y2=y(index(:,2)); y3=y(index(:,3));
invdet=1./(x1.*(y2-y3)-x2.*(y1-y3)+x3.*(y1-y2));

%get alpha and beta
alpha=[invdet.*(y2-y3) invdet.*(y3-y1) invdet.*(y1-y2)];
beta =[invdet.*(x3-x2) invdet.*(x1-x3) invdet.*(x2-x1)];

%get gamma if requested
if nargout==3,
	gamma=zeros(nels,3);
	gamma=[invdet.*(x2.*y3-x3.*y2) invdet.*(y1.*x3-y3.*x1) invdet.*(x1.*y2-x2.*y1)];
	varargout{1}=gamma;
end
