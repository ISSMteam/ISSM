function plot_riftvel(md,options,nlines,ncols,index)
%PLOT_RIFTVEL - plot rift velocity
%
%   Usage:
%      plot_riftvel(md,options,width,i);
%
%   See also: PLOTMODEL

%some checks
if (length(md.initialization.vx)~=md.mesh.numberofvertices | length(md.initialization.vy)~=md.mesh.numberofvertices),
	error('plot_riftvel error message: vx and vy do not have the right size'),
end
if ~isstruct(md.rifts.riftstruct),
	error('plot error message: no rifts available!');
end
options=addfielddefault(options,'scaling',2);

%set as NaN all velocities not on rifts
u=NaN*ones(md.mesh.numberofvertices,1);
v=NaN*ones(md.mesh.numberofvertices,1);
for i=1:size(md.rifts.riftstruct,1),
	penaltypairs=md.rifts.riftstruct(i).penaltypairs(:,[1 2]);
	u(penaltypairs(:))=md.initialization.vx(penaltypairs(:));
	v(penaltypairs(:))=md.initialization.vy(penaltypairs(:));
end

%process data and model
[x y z elements is2d isplanet]=processmesh(md,[],options);
[vel datatype]=processdata(md,[u v],options);
[quivers,palette]=quiver_process(x,y,vel(:,1),vel(:,2),options);

%prepare plot
subplot(nlines,ncols,index); 
hold on

%plot mesh boundaries
for i=1:size(md.mesh.segments,1),
	plot(x(md.mesh.segments(i,1:2)),y(md.mesh.segments(i,1:2)),'k.-');
end

%plot rifts vel
h3=[];
for i=1:quivers.numcolors
	pos=find(quivers.colorind==i);
	hprime=quiver(quivers.x(pos),quivers.y(pos),quivers.u(pos),quivers.v(pos),...
		'Color',palette(i,:),'ShowArrowHead','on','AutoScale','off');
	h3=[h3;hprime];
end

%plot rift velocities
isp1=0;
isp2=0;

%plot mesh boundaries
for i=1:size(md.mesh.segments,1),
	h1=plot(x(md.mesh.segments(i,1:2)),y(md.mesh.segments(i,1:2)),'b-');
end

for i=1:size(md.rifts.riftstruct,1),
	%get nodes on rift
	penaltypairs=md.rifts.riftstruct(i).penaltypairs;

	segments=md.rifts.riftstruct(i).segments;
	for j=1:size(segments,1),
		plot(x(segments(j,1:2)),y(segments(j,1:2)),'b-');
	end

	normal=zeros(2,1);
	for j=1:size(penaltypairs,1),
		normal(1)=penaltypairs(j,5);
		normal(2)=penaltypairs(j,6);

		vx1=md.initialization.vx(penaltypairs(j,1)); vx2=md.initialization.vx(penaltypairs(j,2)); vy1=md.initialization.vy(penaltypairs(j,1)); vy2=md.initialization.vy(penaltypairs(j,2));
		penetration=(vx2-vx1)*normal(1)+(vy2-vy1)*normal(2);
		%if penetration is negative, plot in black, positive, plot in red;: ie: if rift is closing, black, if rift is opening, red.
		if(penetration>0),
			p2=plot(x(penaltypairs(j,1)) ,y(penaltypairs(j,1)),'*'); set(p2,'Color',[140 140 140]/255);
			isp2=1;
		else
			p1=plot(x(penaltypairs(j,1)) ,y(penaltypairs(j,1)),'k*');
			isp1=1;
		end
	end

	%point out the tips
	h2=plot(x(md.rifts.riftstruct(i).tips(1)),y(md.rifts.riftstruct(i).tips(1)),'g*');
	plot(x(md.rifts.riftstruct(i).tips(2)),y(md.rifts.riftstruct(i).tips(2)),'g*');
	segments=md.rifts.riftstruct(i).segments(:,1:2);
end

%legend
if isp1 & isp2
	legend([h1,h2,p1,p2],'mesh boundaries','rift tips',' rifts closing','rifts opening')
elseif isp1
	legend([h1,h2,p1],'mesh boundaries','rift tips',' rifts closing')
elseif isp2
	legend([h1,h2,p2],'mesh boundaries','rift tips','rifts opening')
else
	legend([h1,h2],'mesh boundaries','rift tips')
end
hold off

%apply options
quiver_colorbar(quivers,options);
options=changefieldvalue(options,'colorbar',2);
options=addfielddefault(options,'title','Rift Velocities');
applyoptions(md,[],options);
