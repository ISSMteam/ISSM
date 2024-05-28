function plot_riftnumbering(md,options,nlines,ncols,index)
%PLOT_RIFTNUMBERING - plot rift penetration + numbering of all rift vertices, as well as rift numbers.
%
%   Usage:
%      plot_riftnumbering(md,options,width,i);
%
%   See also: PLOTMODEL

%process data and model
[x y z elements is2d isplanet]=processmesh(md,[],options);
fontsize=getfieldvalue(options,'FontSize',8);

subplot(nlines,ncols,index); 
hold on

%plot mesh boundaries
for i=1:size(md.mesh.segments,1),
	plot(x(md.mesh.segments(i,1:2)),y(md.mesh.segments(i,1:2)),'k.-');
end

isp1=0;
isp2=0;

if isstruct(md.rifts.riftstruct),
	%plot mesh boundaries
	for i=1:size(md.mesh.segments,1),
		h1=plot(x(md.mesh.segments(i,1:2)),y(md.mesh.segments(i,1:2)),'b-');
	end
	for i=1:size(md.rifts.riftstruct,1),
		penaltypairs=md.rifts.riftstruct(i).penaltypairs;

		segments=md.rifts.riftstruct(i).segments;
		for j=1:size(segments,1),
			plot(x(segments(j,1:2)),y(segments(j,1:2)),'b-');
		end

		normal=zeros(2,1);
		for j=1:size(penaltypairs,1),
			normal(1)=penaltypairs(j,5);
			normal(2)=penaltypairs(j,6);

			vx1=md.initialization.vx(penaltypairs(j,1)); 
			vx2=md.initialization.vx(penaltypairs(j,2));
			vy1=md.initialization.vy(penaltypairs(j,1)); 
			vy2=md.initialization.vy(penaltypairs(j,2));
			penetration=(vx2-vx1)*normal(1)+(vy2-vy1)*normal(2);
			%if penetration is negative, plot in black, positive, plot in red;: ie: if rift is closing, black, if rift is opening, red.
			if(penetration>0),
				p2=plot(x(penaltypairs(j,1)) ,y(penaltypairs(j,1)),'ro-','LineWidth',1);
				set(p2,'MarkerSize',3);
				isp2=1;
			else
				p1=plot(x(penaltypairs(j,1)) ,y(penaltypairs(j,1)),'ko-','LineWidth',1);
				set(p1,'MarkerSize',3);
				isp1=1;
			end
		end

		%point out the tips
		h2=plot(x(md.rifts.riftstruct(i).tips(1)),y(md.rifts.riftstruct(i).tips(1)),'g*');
		plot(x(md.rifts.riftstruct(i).tips(2)),y(md.rifts.riftstruct(i).tips(2)),'g*');
	end
	if strcmpi(getfieldvalue(options,'legend','on'),'on'),
		if isp1 & isp2
			l=legend([h1,h2,p1,p2],'mesh boundaries','crack tips','faults','rifts');
		elseif isp1
			l=legend([h1,h2,p1],'mesh boundaries','crack tips','faults');
		elseif isp2
			l=legend([h1,h2,p2],'mesh boundaries','crack tips','rifts');
		else
			l=legend([h1,h2],'mesh boundaries','crack tips');
		end
		set(l,'Location',getfieldvalue(options,'legend_location','NorthEast'));
	end
else
	error('plot error message: no rifts available!');
end

%Now, plot rift vertices numbers.
for i=1:size(md.rifts.riftstruct,1),
	penaltypairs=md.rifts.riftstruct(i).penaltypairs;

	for j=1:size(penaltypairs,1),
		node=penaltypairs(j,1);
		t=text(x(node),y(node),[num2str(i) '.' num2str(j)]);
		set(t,'FontSize',fontsize);
	end
end

%apply options
options=addfielddefault(options,'title','Rift/Fault location');
options=addfielddefault(options,'colorbar',0);
applyoptions(md,[],options);
