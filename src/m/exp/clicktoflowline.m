function clicktoflowline(index,x,y,u,v,filename)
%CLICKTOFLOWLINE - create a flowline ARGUS file
%
%   create a flowline contour file (named 'filename') by clicking
%   on a velocity field once (velocity must be plotted first)
%
%   Usage: 
%      clicktoflowline(index,x,y,u,v,x0,y0,filename)
%
%   Example: 
%      clicktoflowline(md.mesh.elements,md.mesh.x,md.mesh.y,md.inversion.vx_obs,md.inversion.vy_obs,'flowline.exp')

%Get click position
[x0,y0]=exp_ginput(1,pairoptions());

%Get flowline
line=flowlines(index,x,y,u,v,x0,y0);

%plot
hold on
plot(line.x,line.y,'r-');

%Write argus file
expwrite(line,filename);
