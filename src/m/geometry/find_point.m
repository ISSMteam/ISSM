function f=find_point(tabx,taby,pointx,pointy)
%FIND_POINT - find closest point
%
%   find which point of the list (tabx,taby) is the closest to (pointx,pointy)
%
%   Usage:
%      f=find_point(tabx,taby,pointx,pointy)

%Compute distance between point and cloud of points
distance=sqrt((tabx-pointx).^2+(taby-pointy).^2);

%find index of the minimum distance and return the first one only
%f=find(distance==min(min(distance)),1);
[~,f]=min(distance(:));
