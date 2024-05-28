%
%  create an edge perimeter list for the elements in the model.
%
%  [edgeper,elemper,iloop]=edgeperimeter(elem,nodecon,edgeadj)
%
%  where the required input is:
%    elem          (numeric, element connectivity array (elems x nodes))
%    nodecon       (numeric, node connectivity array (nodes x elems+1))
%
%  and the required output is:
%    edgeper       (numeric, edge perimeter list (edgeper x 2))
%    elemper       (numeric, element perimeter list (edgeper x 1))
%    iloop         (numeric, index for each loop (nloop))
%
%  the optional input is:
%    edgeadj       (numeric, edge adjacency array (elems x edges))
%
function [edgeper,elemper,iloop]=edgeperimeter(elem,nodecon,edgeadj)

if ~nargin
    help edgeperimeter
    return
end

%%  create the edge adjacency array

if ~exist('edgeadj','var') || isempty(edgeadj)
    edgeadj=edgeadjacency(elem,nodecon);
end

%%  create the unshared edge list

[icol,irow]=find(edgeadj'==0);
edgeuns=zeros(length(irow),2);
elemuns=zeros(length(irow),1);

%  loop over the edges

for i=1:length(irow)
    edgeuns(i,1)=elem(irow(i),icol(i));
    edgeuns(i,2)=elem(irow(i),mod(icol(i),size(elem,2))+1);
    elemuns(i)=irow(i);
end

%%  create the edge perimeter list

edgeper=zeros(size(edgeuns));
elemper=zeros(size(elemuns));
iloop=[];
ipt=0;

%  find the beginning of a loop

while ~isempty(find(edgeuns,1))
    ipt=ipt+1;
    iloop(end+1)=ipt;
    [irow,icol]=find(edgeuns,1);
    edgeper(ipt,:)=edgeuns(irow,:);
    elemper(ipt)  =elemuns(irow);
    edgeuns(irow,:)=[0 0];
    elemuns(irow)  =0;
    [irow,icol]=find(edgeuns==edgeper(ipt,2),1);

%  continue following the loop

    while ~isempty(irow)
        ipt=ipt+1;
        if (icol == 1)
            edgeper(ipt,:)=edgeuns(irow,:);
        else
            edgeper(ipt,1)=edgeuns(irow,2);
            edgeper(ipt,2)=edgeuns(irow,1);
        end
        elemper(ipt)  =elemuns(irow);
        edgeuns(irow,:)=[0 0];
        elemuns(irow)  =0;
        [irow,icol]=find(edgeuns==edgeper(ipt,2),1);
    end

%  check to see if loop is closed

    if (edgeper(iloop(end),1) ~= edgeper(ipt,2))
        warning('Loop %d is not closed.\n',length(loop));
    end
end

end
