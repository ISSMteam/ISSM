function plot_sarpwr(md,options,width,i)
%PLOT_SARPWR - plot radar image
%
%   Usage:
%      plot_sarpwr(md,options,width,i);
%
%   See also: PLOTMODEL

%plot mesh sarpwr
subplot(width,width,i); 

%units
if exist(options,'unit'),
	unit=getfieldvalue(options,'unit');
	md.mesh.x=md.mesh.x*unit;
	md.mesh.y=md.mesh.y*unit;
	md.mesh.z=md.mesh.z*unit;
end

imagesc(md.radaroverlay.x,md.radaroverlay.y,double(md.radaroverlay.pwr)),set(gca,'YDir','normal');colormap(gray);

%apply options
options=addfielddefault(options,'colorbar',0);
options=changefieldvalue(options,'colormap','gray');

applyoptions(md,[],options);
