function plot_boundaries(md,options,width,i)
%PLOT_BOUNDARIES - plot mesh boundaries
%
%   Usage:
%      plot_boundaries(md,options,width,i);
%
%   See also: PLOTMODEL

subplot(width,width,i); 

%process data and model
if getfieldvalue(options,'layer',0)
	options=removefield(options,'layer',1);
end
[x y z elements is2d isplanet]=processmesh(md,[],options);

for i=1:size(md.mesh.segments,1),
	plot(x(md.mesh.segments(i,1:2)),y(md.mesh.segments(i,1:2)),'k.-');hold on;
end

%plot rifts if present: 
if isstruct(md.rifts.riftstruct),
	for i=1:size(md.rifts.riftstruct,1),
		segments=md.rifts.riftstruct(i).segments;
		for j=1:size(segments,1),
			plot(x(segments(j,1:2)),y(segments(j,1:2)),'r.-');
		end
		text(x(segments(floor(size(segments,1)/4),1)),y(segments(floor(size(segments,1)/4),1)),['Rift #' num2str(i)]);
		%point out the tips
		plot(x(md.rifts.riftstruct(i).tips(1)),y(md.rifts.riftstruct(i).tips(1)),'b*');
		plot(x(md.rifts.riftstruct(i).tips(2)),y(md.rifts.riftstruct(i).tips(2)),'b*');
	end
end

%apply options
options=addfielddefault(options,'title','Mesh boundaries');
options=addfielddefault(options,'colorbar',0);
options=addfielddefault(options,'view',2);
applyoptions(md,[],options);
