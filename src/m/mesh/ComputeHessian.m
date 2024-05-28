function hessian=ComputeHessian(index,x,y,field,type)
%COMPUTEHESSIAN - compute hessian matrix from a field
%
%   Compute the hessian matrix of a given field
%   return the three components Hxx Hxy Hyy
%   for each element or each node
%
%   Usage:
%      hessian=ComputeHessian(index,x,y,field,type)
%
%   Example:
%      hessian=ComputeHessian(md.mesh.elements,md.mesh.x,md.mesh.y,md.inversion.vel_obs,'node')

%some variables
numberofnodes=length(x);
numberofelements=size(index,1);

%some checks
if length(field)~=numberofnodes & length(field)~=numberofelements,
	error('ComputeHessian error message: the given field size not supported yet');
end
if ~strcmpi(type,'node') & ~strcmpi(type,'element'),
	error('ComputeHessian error message: only ''node'' or ''element'' type supported yet');
end

%initialization
line=index(:);
linesize=3*numberofelements;

%get areas and nodal functions coefficients N(x,y)=alpha x + beta y + gamma 
[alpha beta]=GetNodalFunctionsCoeff(index,x,y);
areas=GetAreas(index,x,y);

%compute weights that hold the volume of all the element holding the node i
weights=sparse(line,ones(linesize,1),repmat(areas,3,1),numberofnodes,1);

%compute field on nodes if on elements
if length(field)==numberofelements,
	field=sparse(line,ones(linesize,1),repmat(areas.*field,3,1),numberofnodes,1)./weights ;
end

%Compute gradient for each element
grad_elx=sum(field(index).*alpha,2); 
grad_ely=sum(field(index).*beta,2);

%Compute gradient for each node (average of the elements around)
gradx=sparse(line,ones(linesize,1),repmat(areas.*grad_elx,3,1),numberofnodes,1);
grady=sparse(line,ones(linesize,1),repmat(areas.*grad_ely,3,1),numberofnodes,1);
gradx=gradx./weights;
grady=grady./weights;

%Compute hessian for each element
hessian=[sum(gradx(index).*alpha,2) sum(grady(index).*alpha,2) sum(grady(index).*beta,2)];

if strcmpi(type,'node')
	%Compute Hessian on the nodes (average of the elements around)
	hessian=[sparse(line,ones(linesize,1),repmat(areas.*hessian(:,1),3,1),numberofnodes,1)./weights ...
		sparse(line,ones(linesize,1),repmat(areas.*hessian(:,2),3,1),numberofnodes,1)./weights ...
		sparse(line,ones(linesize,1),repmat(areas.*hessian(:,3),3,1),numberofnodes,1)./weights ];
end
