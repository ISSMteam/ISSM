function plot_riftrelvel(md,options,nlines,ncols,index)
%PLOT_RIFTRELVEL - plot rift relative velocities
%
%   Usage:
%      plot_riftrelvel(md,options,nlines,ncols,i);
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

%markersize: 
markersize=getfieldvalue(options,'markersize',1);

%recover vx and vy:
vx=getfieldvalue(options,'riftrelvel_vx',md.initialization.vx);
vy=getfieldvalue(options,'riftrelvel_vy',md.initialization.vy);

%set as NaN all velocities not on rifts
u=NaN*ones(md.mesh.numberofvertices,1);
v=NaN*ones(md.mesh.numberofvertices,1);
for i=1:size(md.rifts.riftstruct,1),
	penaltypairs=md.rifts.riftstruct(i).penaltypairs(:,[1 2]);
	u(md.rifts.riftstruct(i).penaltypairs(:,1))=vx(penaltypairs(:,1))-vx(penaltypairs(:,2));
	v(md.rifts.riftstruct(i).penaltypairs(:,1))=vy(penaltypairs(:,1))-vy(penaltypairs(:,2));
end

%process data and model
[x y z elements is2d isplanet]=processmesh(md,[],options);
[vel datatype]=processdata(md,[u v],options);
[quivers,palette]=quiver_process(x,y,vel(:,1),vel(:,2),options);

%prepare plot
subplot(nlines,ncols,index); 
hold on

%plot rifts vel
h3=[];
for i=1:quivers.numcolors
	pos=find(quivers.colorind==i);
	hprime=quiver(quivers.x(pos),quivers.y(pos),quivers.u(pos),quivers.v(pos),...
		'Color',palette(i,:),'ShowArrowHead','on','AutoScale','off');
	hprime=quiver(quivers.x(pos),quivers.y(pos),-quivers.u(pos),-quivers.v(pos),...
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
		plot(x(segments(j,1:2)),y(segments(j,1:2)),'k-');
	end

	normal=zeros(2,1);
	for j=1:size(penaltypairs,1),
		normal(1)=penaltypairs(j,5);
		normal(2)=penaltypairs(j,6);

		vx1=vx(penaltypairs(j,1)); vx2=vx(penaltypairs(j,2)); vy1=vy(penaltypairs(j,1)); vy2=vy(penaltypairs(j,2));
		penetration=(vx2-vx1)*normal(1)+(vy2-vy1)*normal(2);
		%if penetration is negative, plot in black, positive, plot in red;: ie: if rift is closing, black, if rift is opening, red.
		if(penetration>0),
			p2=plot(x(penaltypairs(j,1)) ,y(penaltypairs(j,1)),'.','MarkerSize',markersize); set(p2,'Color',[140 140 140]/255);
			isp2=1;
		else
			p1=plot(x(penaltypairs(j,1)) ,y(penaltypairs(j,1)),'k.','MarkerSize',markersize);
			isp1=1;
		end
	end

	%point out the tips
	h2=plot(x(md.rifts.riftstruct(i).tips(1)),y(md.rifts.riftstruct(i).tips(1)),'g.','MarkerSize',markersize);
	plot(x(md.rifts.riftstruct(i).tips(2)),y(md.rifts.riftstruct(i).tips(2)),'g.','MarkerSize',markersize);
	segments=md.rifts.riftstruct(i).segments(:,1:2);
end

faulttitle=getfieldvalue(options','faulttitle','faults');
rifttitle=getfieldvalue(options','rifttitle','rifts');
%legend
if strcmpi(getfieldvalue(options,'legend','on'),'on'),
	if isp1 & isp2
		l=legend([h1,h2,p1,p2],'mesh boundaries','crack tips',faulttitle,rifttitle);
	elseif isp1
		l=legend([h1,h2,p1],'mesh boundaries','crack tips',faulttitle);
	elseif isp2
		l=legend([h1,h2,p2],'mesh boundaries','crack tips',rifttitle);
	else
		l=legend([h1,h2],'mesh boundaries','crack tips');
	end
set(l,'Location',getfieldvalue(options,'legend_location','NorthEast'));
end
hold off

%apply options
quiver_colorbar(quivers,options);
options=changefieldvalue(options,'colorbar',2);
options=addfielddefault(options,'title','Rift/Fault Relative Velocity');
applyoptions(md,[],options);
