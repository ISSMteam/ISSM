function plot_section(md,data,options,nlines,ncols,i)
%PLOT_SECTION - plot a given field on a section
%
%   Usage:
%      plot_section(md,data,options,nlines,ncols,i)
%
%   See also: PLOTMODEL

%How many subplots?
if exist(options,'showsection')

	%Compute the indexes of the 2 plots (one for the sectionvalue and one for showsection
	upperplots=floor((i-1)/ncols);
	if upperplots==0, leftplots=i-1; else leftplots=i-ncols*upperplots-1; end
	index1=4*ncols*upperplots+2*leftplots+1;
	index2=index1+1;
	ncols=2*ncols;
else
	index1=i;
end

%process model
[x_m y_m z_m elements_m is2d isplanet]=processmesh(md,[],options);

%Get number of curves and generate random colors
numcurves=size(data,2);
colorm=getfieldvalue(options,'colormap','lines');
color=eval([ colorm '(numcurves);']);
options=removefield(options,'colormap',0); %back to default colormap

%replug x and y onto model so that SectionValue treats the problem correctly
md3d=md;
if exist(options,'layer')
	md.mesh.x=md.mesh.x2d; md.mesh.y=md.mesh.y2d; md.mesh.elements=md.mesh.elements2d;
	md.mesh=mesh2d(md.mesh);
end

%read contours: 
profiles=expread(getfieldvalue(options,'sectionvalue'));
numprofiles=length(profiles);

%Loop over number of profiles: 
for profile_i=1:numprofiles,
	profile=profiles(profile_i);

	%Loop over number of curves
	for i=1:numcurves,

		[datai datatype]=processdata(md3d,data(:,i),options);

		%resolution
		if exist(options,'resolution'),
			resolution=getfieldvalue(options,'resolution');
		else %Default resolution
			if is2d,
				resolution=[1000 1];
			else
				resolution=[1000 10*md.mesh.numberoflayers];
			end
			disp(['plot_section warning: no resolution specified, use default resolution: [horizontal_resolution vertical_resolution]=[' num2str(resolution)  ']']);
		end

		%Compute section value
		[elements,x,y,z,s,data_s]=SectionValues(md,datai,profile,resolution);

		if getfieldvalue(options,'sectionmean',0)==1,
			disp(['Mean value of data along section: ' num2str(mean(data_s))])
			disp(['Median value of data along section: ' num2str(median(data_s))])
			disp(['Standard deviation of data along section: ' num2str(std(data_s))])
		end

		%units
		if exist(options,'unit'),
			unit=getfieldvalue(options,'unit');
			x=x*unit;
			y=y*unit;
			z=z*unit;
			s=s*unit;
		end

		%2D
		if is2d,
%		%plot section value
%		hold on;
%		subplot(nlines,ncols,index1)
%		%subplot(1,3,[2 3])
%		plot(s,data_s,'color',color(i,:),'LineWidth',getfieldvalue(options,'linewidth',1))
%		%3D
%	else
%		%plot section value
%		%if user requested view2: 2d plot with curvilinear coordinate
%		if (getfieldvalue(options,'view',3)==2 )

			%Show Section if requested by user
			if exist(options,'showsection')

				%compute number of labels
				numlabels=min(getfieldvalue(options,'showsection'),length(s));
				shift=fix(length(s)/numlabels);

				%plot labels on current graph
				hold on
				text(s(1),data_s(1),'1','backgroundcolor',[0.8 0.9 0.8])
				for i=2:numlabels-1
					text(s(1+(i-1)*shift),data_s(1+(i-1)*shift),num2str(i),'backgroundcolor',[0.8 0.9 0.8])
				end
				text(s(end),data_s(end),'end','backgroundcolor',[0.8 0.9 0.8])

				%plot section only with labels
				subplot(nlines,ncols,index2)
				plot_unit(x_m,y_m,z_m,elements_m,data(:,i),is2d,isplanet,datatype,options)
				hold on
				text(x(1),y(1),'1','backgroundcolor',[0.8 0.9 0.8])
				for i=2:numlabels-1
					text(x(1+(i-1)*shift),y(1+(i-1)*shift),num2str(i),'backgroundcolor',[0.8 0.9 0.8])
				end
				text(x(end),y(end),'end','backgroundcolor',[0.8 0.9 0.8])
				plot(x,y,'-r')
				axis([min(md.mesh.x)-1 max(md.mesh.x)+1 min(md.mesh.y)-1 max(md.mesh.y)+1])
				view(2)
			end

			%plot section value
			if(i==1), subplot(nlines,ncols,index1); end
			plot(s,data_s,'color',color(i,:),'LineWidth',getfieldvalue(options,'linewidth',1))
			hold on

			%3D
		else
			%plot section value
			%if user requested view2: 2d plot with curvilinear coordinate
			if (getfieldvalue(options,'view',3)==2 )

				%Show Section if requested by user
				if exist(options,'showsection')

					%compute number of labels
					numlabels=min(getfieldvalue(options,'showsection'),length(s));
					shift=fix(length(s)/numlabels);

					%plot labels on current graph
					hold on
					text(s(1),z(1),'1','backgroundcolor',[0.8 0.9 0.8])
					for i=2:numlabels-1
						text(s(1+(i-1)*shift),z(1+(i-1)*shift),num2str(i),'backgroundcolor',[0.8 0.9 0.8])
					end
					text(s(end),z(end),'end','backgroundcolor',[0.8 0.9 0.8])

					%plot section only with labels
					subplot(nlines,ncols,index2)
					plot_unit(x_m,y_m,z_m,elements_m,data(:,i),is2d,datatype,options)
					hold on
					text(x(1),y(1),'1','backgroundcolor',[0.8 0.9 0.8])
					for i=2:numlabels-1
						text(x(1+(i-1)*shift),y(1+(i-1)*shift),num2str(i),'backgroundcolor',[0.8 0.9 0.8])
					end
					text(x(end),y(end),'end','backgroundcolor',[0.8 0.9 0.8])
					plot(x,y,'-r')
					axis([min(md.mesh.x)-1 max(md.mesh.x)+1 min(md.mesh.y)-1 max(md.mesh.y)+1])
					view(2)
				end

				subplot(nlines,ncols,index1)
				A=elements(:,1); B=elements(:,2); C=elements(:,3);  D=elements(:,4); 
				patch( 'Faces', [A B C D], 'Vertices', [s z zeros(length(s),1)],'FaceVertexCData',data_s,'FaceColor','interp','EdgeColor','none');

			else

				%Show Section if requested by user
				if exist(options,'showsection')

					%compute number of labels
					numlabels=min(getfieldvalue(options,'showsection'),length(s));
					shift=fix(length(x)/numlabels);

					%plot labels on current graph
					hold on
					text(x(1),y(1),z(1),'1','backgroundcolor',[0.8 0.9 0.8])
					for i=2:numlabels-1
						text(x(1+(i-1)*shift),y(1+(i-1)*shift),z(1+(i-1)*shift),num2str(i),'backgroundcolor',[0.8 0.9 0.8])
					end
					text(x(end),y(end),z(end),'end','backgroundcolor',[0.8 0.9 0.8])

					%plot section only with labels
					subplot(nlines,ncols,index2)
					plot_unit(x_m,y_m,z_m,elements_m,data,is2d,datatype,options)
					hold on
					text(x(1),y(1),'1','backgroundcolor',[0.8 0.9 0.8])
					for i=2:numlabels-1
						text(x(1+(i-1)*shift),y(1+(i-1)*shift),num2str(i),'backgroundcolor',[0.8 0.9 0.8])
					end
					text(x(end),y(end),'end','backgroundcolor',[0.8 0.9 0.8])
					plot(x,y,'-r')
					axis([min(md.mesh.x)-1 max(md.mesh.x)+1 min(md.mesh.y)-1 max(md.mesh.y)+1])
					view(2)
				end

				subplot(nlines,ncols,index1)
				A=elements(:,1); B=elements(:,2); C=elements(:,3);  D=elements(:,4); 
				patch( 'Faces', [A B C D], 'Vertices', [x y z],'FaceVertexCData',data_s,'FaceColor','interp','EdgeColor','none');
				view(3)
			end
		end
	end
end

%apply options
options=addfielddefault(options,'title','Section value');
if dimension(md.mesh)==2
	options=addfielddefault(options,'colorbar',0);
end
if ((dimension(md.mesh)==2) | getfieldvalue(options,'view',2)==2 )
	options=addfielddefault(options,'xlabel','Curvilinear coordinate');
	options=addfielddefault(options,'axis','auto');
end
if (dimension(md.mesh)==3 & getfieldvalue(options,'view',2)==2 )
	options=addfielddefault(options,'ylabel','z');
	options=addfielddefault(options,'axis','auto');
end
applyoptions(md3d,[],options);
