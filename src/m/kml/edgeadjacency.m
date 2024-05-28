%
%  create an edge adjacency table for the elements in the model.
%
%  [edgeadj]=edgeadjacency(elem,nodecon)
%
%  where the required input is:
%    elem          (numeric, element connectivity array (elems x nodes))
%    nodecon       (numeric, node connectivity array (nodes x elems+1))
%
%  and the required output is:
%    edgeadj       (numeric, edge adjacency array (elems x edges))
%
function [edgeadj]=edgeadjacency(elem,nodecon)

if ~nargin
    help edgeadjacency
    return
end

%%  create the edge adjacency array

edgeadj=zeros(size(elem));

%  loop over the elements

for i=1:size(elem,1)

%  loop over the edges for each element (trias only for now)

    for j=1:size(elem,2)
        inode1=elem(i,j);
        inode2=elem(i,mod(j,size(elem,2))+1);

%  loop over the elements containing the first node of the edge to see
%  if they contain the second node of the edge

        for k=1:nodecon(inode1,end)
            if (nodecon(inode1,k) ~= i) && ...
               ~isempty(find(elem(nodecon(inode1,k),:)==inode2,1))
                edgeadj(i,j)=nodecon(inode1,k);
                break;
            end
        end
    end
end

end
