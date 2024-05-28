function [gradx, grady]=computeGrad(index,x,y,field)
%COMPUTEGRAD - compute the gradient from a field

%some variables
numberofnodes=length(x);
numberofelements=size(index,1);
numberoftime = size(field, 2);

%some checks
if (length(field)~=numberofnodes) && (length(field)~=numberofelements)
    error('ComputeHessian error message: the given field size not supported yet');
end
%initialization
line=index(:);
linesize=3*numberofelements;

%get areas and nodal functions coefficients N(x,y)=alpha x + beta y + gamma
[alpha, beta]=GetNodalFunctionsCoeff(index,x,y);
areas=GetAreas(index,x,y);

%compute weights that hold the volume of all the element holding the node i
weights=sparse(line,ones(linesize,1),repmat(areas,3,1),numberofnodes,1);

%compute field on nodes if on elements
if length(field)==numberofelements
    field=sparse(line,ones(linesize,1),repmat(areas.*field,3,1),numberofnodes,1)./weights ;
end

%Compute gradient for each element
if numberoftime == 1
    grad_elx=sum(field(index).*alpha,2);
    grad_ely=sum(field(index).*beta,2);
else
    grad_elx = zeros(numberofelements,numberoftime);
    grad_ely = zeros(numberofelements,numberoftime);    
    for i = 1:3
        grad_elx = grad_elx + field(index(:,i), :).*alpha(:,i);
        grad_ely = grad_ely + field(index(:,i), :).*beta(:,i);
    end
end
%Compute gradient for each node (average of the elements around)
gradx=sparse(repmat(line,1,numberoftime),cumsum(ones(linesize,numberoftime),2),repmat(areas.*grad_elx,3,1),numberofnodes,numberoftime);
grady=sparse(repmat(line,1,numberoftime),cumsum(ones(linesize,numberoftime),2),repmat(areas.*grad_ely,3,1),numberofnodes,numberoftime);
gradx=gradx./weights;
grady=grady./weights;
