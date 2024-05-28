function plot_riftfraction(md,options,nlines,ncols,index)
%PLOT_RIFTFRACTION - plot rift fractions
%
%   Usage:
%      plot_riftfraction(md,options,width,i);
%
%   See also: PLOTMODEL

%check that there is something in riftproperties
if isnan(md.rifts.riftstruct.riftproperties),
	error('plot_riftfraction error message: field riftproperies is empty, run the model first')
end

%process data and model
[x y z elements is2d isplanet]=processmesh(md,[],options);

subplot(nlines,ncols,index); 
hold on

%plot mesh boundaries
for i=1:size(md.mesh.segments,1),
	h1=plot(x(md.mesh.segments(i,1:2)),y(md.mesh.segments(i,1:2)),'k.-');
end

%first, build a vector of fractions, over all nodes. 
fractions=zeros(md.mesh.numberofvertices,1);

%complete the tips.
for i=1:length(md.rifts.riftstruct), 
	tips=md.rifts.riftstruct(i).tips;
	fractions(tips)=1;
end

hold on;
for i=1:length(md.rifts.riftstruct), 
	segments=md.rifts.riftstruct(i).segments(:,1:2)';
	xc=x(segments(:));
	yc=y(segments(:));
	zc=fractions(segments(:));
	h2=patch('Xdata',xc,'Ydata',yc,'Zdata',zc,'Cdata',zc,'facecolor','none','edgecolor','flat');
end
legend([h1,h2],'mesh boundaries','rifts')
hold off

%apply options
options=addfielddefault(options,'title','Rift ice/water fraction ???????'); %Eric, could you enter a better title?
options=addfielddefault(options,'colorbar',0);
applyoptions(md,[],options);
