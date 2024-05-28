function avg=zoneaverage(md,field,parameter,width,distweight,areaweight,iter)
%ZONEAVERAGE - averages the field over the elements within the distance specified by parameter .* width. 
%               Average is weighted by the element area using model's triangulation.
%               Width determines distance multiplier of parameter when calculating a node's distance for inclusion in an average.
%               Distweight and Areaweight control how the node distance and area affect the weighted average.
%               Iteration is the number of iterations to run the zone averaging.
%
%   Usage:
%      avg=zone_average(md,field,parameter,width,distweight,areaweight)

%TODO: HANDLE NARGINS/PAIROPTION STYLE ARGS FOR WIDTH, DISTWEIGHT, AREAWEIGHT, ITER

if (dimension(md.mesh)==2),
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

%node
if numel(field)==numberofnodes,
    areas=GetAreas(md.mesh.elements,md.mesh.x,md.mesh.y);

    %calculate distances between nodes and reshape into addressable matrix
    nodes=[md.mesh.x, md.mesh.y];
    dist=pdist([nodes;nodes]);
    dist=squareform(dist);
    dist=dist(1:numberofnodes,numberofnodes+1:end);

    %calculate max distance to average each vertex using parameter * width and mask distances
    maxdist = parameter * width;
    flags=dist<=maxdist;
    
    %calculate average area per node
    nodeareas=zeros(numberofnodes,1);
    elementspernode=zeros(numberofnodes,1);
    for i=1:numberofelements,
        for j=1:3,
            node=index(i,j);
            nodeareas(node)=nodeareas(node)+areas(i);
            elementspernode(node)=elementspernode(node)+1;
        end
    end
    nodeareas=nodeareas./elementspernode;
    
    %get distance/area weighted average of field within specified parameter width
    avg=zeros(numberofnodes,1);
    totalweight=distweight+areaweight;
    distweight=distweight/totalweight;
    areaweight=areaweight/totalweight;
    for iter=1:iter,
%     for iter=1:1,
        for i=1:numberofnodes,
            %distance weighting
            distsum = sum(dist(i,flags(i,:))) + eps;
            distmax = max(dist(i,flags(i,:)));
            distweights = ones(sum(flags(i,:)),1);
            distweights = distweights + 1 - distweight*(dist(i,(flags(i,:)))'/distsum); %inverse linear distance weight
            distweights = distweights / sum(distweights);

            %area weighting
            areasum = sum(nodeareas(flags(i,:))) + eps;
            areaweights = areaweight*nodeareas(flags(i,:))/areasum;
            weights=distweights+areaweights;
            avg(i)=nansum(weights.*field(flags(i,:)));
            %avg(i)=nanmean(field(flags(i,:)),1);

    %         areasum = sum(nodeareas(flags));
    %         weights = nodeareas(flags)/areasum;
    %         avg=nanmean(weights.*field(flags),2);
        end
        field=avg;
    end
end
