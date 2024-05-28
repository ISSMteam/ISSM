function expcircle(filename,x0,y0,radius,numberofnodes)
%EXPCIRCLE - create a circular contour corresponding to given parameters
%
%   Creates a closed argus contour centered on x,y of radius size.
%   The contour is made of numberofnodes
%
%   Usage:
%      expcircle(filename,x0,y0,radius,numberofnodes)
%
%   See also EXPMASTER, EXPDOC

%Calculate the cartesians coordinates of the points
x_list=ones(numberofnodes+1,1);
y_list=ones(numberofnodes+1,1);

theta=(0:2*pi/numberofnodes:2*pi*(1-1/numberofnodes))';
theta=[theta;0];

x_list=radius*x_list.*cos(theta);
y_list=radius*y_list.*sin(theta);

%offset x_list and y_list by x0 and y0:
x_list=x_list+x0;
y_list=y_list+y0;

contour.x=x_list;
contour.y=y_list;
contour.density=1;
contour.name='circle';

expwrite(contour,filename);
