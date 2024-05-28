function h=pcolorcen(x,y,z)

% h = pcolorcen(x,y,z);
% h = pcolorcen(z);
%
% like pcolor.m, except it fixes the problem wherein the last row and column of z
% don't get displayed and values are plotted above and to the right of the corresponding
% values in x and y. "cen" is for centered.
%
% in other words: use pcolor when you're specifying the corners of your grid cells, and
% pcolorcen when you're specifying the centers.
%
% x and y can be vectors; pcolorcen will call meshgrid to save you the trouble.
%
% sets shading flat for good measure.
%
% neil banas jul 2002
% (neil@ocean.washington.edu)

if nargin==1
	z = x;
	x = 1:size(z,2);
	y = 1:size(z,1);
end
if ~isempty(find(size(x)==1))
	[x,y] = meshgrid(x,y);
end

% add a border to the x values
xnew = [   x(:,1)-(x(:,2)-x(:,1))            x      x(:,end)+(x(:,end)-x(:,end-1))      ];
xnew = [xnew(1,:)-(xnew(2,:)-xnew(1,:));  xnew;  xnew(end,:)+(xnew(end,:)-xnew(end-1,:))];
% interpolate to add the grid centers
xnew = interp2(xnew,1);
% keep only the grid centers
xnew = xnew(2:2:end,2:2:end);

% repeat for y
ynew = [   y(:,1)-(y(:,2)-y(:,1))            y      y(:,end)+(y(:,end)-y(:,end-1))      ];
ynew = [ynew(1,:)-(ynew(2,:)-ynew(1,:));  ynew;  ynew(end,:)+(ynew(end,:)-ynew(end-1,:))];
ynew = interp2(ynew,1);
ynew = ynew(2:2:end,2:2:end);

% add a border of nans to z
znew = [[z nan.*ones(size(z,1),1)]; nan.*ones(1,size(z,2)+1)];

% now pcolor
h = pcolor(xnew,ynew,znew);
shading flat;