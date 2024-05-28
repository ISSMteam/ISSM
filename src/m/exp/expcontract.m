function normal_node=expcontract(newfile,oldfile,distance)
%EXPCONTRACT - contract or expand a profile, according to the normal.
% 
%   Usage:
%      expcontract(newfile,oldfile,distance)
%
%   See also EXPMASTER, EXPDOC

contour=expread(oldfile);
num=numel(contour.x);

normal=zeros(num-1,2);
normal_node=zeros(num-1,2);

for i=1:num-1,
	normal(i,:)=[ contour.y(i)-contour.y(i+1) contour.x(i+1)-contour.x(i)];
	normal(i,:)=normal(i,:)/sqrt(normal(i,1)^2+normal(i,2)^2);
end

normal_node(2:end,:)=[normal(1:end-1,:)+normal(2:end,:)];
normal_node(1,:)=normal(1,:)+normal(end,:);

normal_node_norm=sqrt(normal_node(:,1).^2+normal_node(:,2).^2);
normal_node(:,1)=normal_node(:,1)./normal_node_norm;
normal_node(:,2)=normal_node(:,2)./normal_node_norm;

contour.x(1:end-1)=contour.x(1:end-1)+distance*normal_node(:,1);
contour.y(1:end-1)=contour.y(1:end-1)+distance*normal_node(:,2);

contour.x(end)=contour.x(1);
contour.y(end)=contour.y(1);

expwrite(contour,newfile);
