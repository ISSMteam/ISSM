function plot_parthist(md,options,nlines,ncols,i)
%PLOT_PARTHIST - plot partitioning histogram
%
%   Usage:
%      plot_parthist(md,options,nlines,ncols,i);
%
%   See also: PLOTMODEL

%process data and model
[x y z elements is2d isplanet]=processmesh(md,[],options);

%plot mesh
subplot(nlines,ncols,i); 

imin=min(md.qmu.partition);
imax=max(md.qmu.partition);

part=zeros(imax-imin+1,2);

for i=imin:imax
    ind=find(md.qmu.partition == i);
    part(i-imin+1,1)=length(ind);
	part(i-imin+1,2)=sum(md.vertex_weight(ind));
end

subplot(2,1,1)
bar(imin:imax,part(:,1));
xlim([imin-0.5 imax+0.5])
title('Number of Nodes in Each Partition')

subplot(2,1,2)
bar(imin:imax,part(:,2));
xlim([imin-0.5 imax+0.5])
title('Total Weight in Each Partition')
