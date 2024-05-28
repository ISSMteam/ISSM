%
%  create a node connectivity table for the elements in the model.
%
%  [nodecon]=edgeadjacency(elem,nnodes,mxepg)
%
%  where the required input is:
%    elem          (numeric, element connectivity array (elems x nodes))
%
%  and the required output is:
%    nodecon       (numeric, node connectivity array (nnodes x mxepg+1))
%
%  the optional input is:
%    nnodes        (numeric, number of nodes)
%    mxepg         (numeric, max elements per node)
%
function [nodecon]=kmlnodeconnectivity(elem,nnodes,mxepg)

if ~nargin
    help kmlnodeconnectivity
    return
end

if ~exist('nnodes','var') || isempty(nnodes)
    nnodes=max(max(elem));
end
if ~exist('mxepg','var') || isempty(mxepg)
    mxepg=25;
end

%%  create the node connectivity array

nodecon=zeros(nnodes,mxepg+1);

%  loop over the elements

for i=1:size(elem,1)

%  loop over the nodes for each element

    for j=1:size(elem,2)
        if elem(i,j)
            nodecon(elem(i,j),nodecon(elem(i,j),end)+1)=i;
            nodecon(elem(i,j),end)=nodecon(elem(i,j),end)+1;
        end
    end
end

%%  sort the node connectivity array

%  loop over the nodes

for i=1:size(nodecon,1)
    if (nodecon(i,end) > 1)
        nodecon(i,1:nodecon(i,end))=sort(nodecon(i,1:nodecon(i,end)));
    end
end

end
