%Test Name: SquareShelfTranSawTooth3d
printingflag = false;

md=triangle(model(),'../Exp/Square.exp',350000.);
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareShelf.par');
md=setflowequation(md,'SSA','all');
md=extrude(md,3,1.);
md.cluster=generic('name',oshostname(),'np',1);
md.transient.isthermal=0;

md.timestepping.time_step=1.;
md.settings.output_frequency=1;
md.timestepping.final_time=2000.;

%Solve for thinning rate -> -1 * surface mass balance
smb= 2.*ones(md.mesh.numberofvertices,1);
md.smb.mass_balance= smb;
md.basalforcings.groundedice_melting_rate= smb;

md=solve(md,'Masstransport');

for i=1:10
	 md=solve(md,'Masstransport');
	 md.smb.mass_balance= md.smb.mass_balance - ((md.results.MasstransportSolution.Thickness)-md.geometry.thickness);
end

%Set up transient
smb = md.smb.mass_balance;

tooth= [ [ones(400,1)*(smb') - 10.]' [ones(400,1)*(smb')]' ];
smb=[ [ones(399,1)*(smb')]' smb  tooth tooth];

md.smb.mass_balance= smb;
md.smb.mass_balance(end+1,:)=[1.:2000.];

md=solve(md,'Transient');

%Fields and tolerances to track changes
field_names={'Vx1','Vy1','Vz1','Vel1','Pressure1','Bed1','Surface1','Thickness1','SmbMassBalance1', ...
	'Vx2','Vy2','Vz2','Vel2','Pressure2','Bed2','Surface2','Thickness2','SmbMassBalance2', ...
	'Vx3','Vy3','Vz3','Vel3','Pressure3','Bed3','Surface3','Thickness3','SmbMassBalance3', ...
	'Vx4','Vy4','Vz4','Vel4','Pressure4','Bed4','Surface4','Thickness4','SmbMassBalance4', ...
	'Vx5','Vy5','Vz5','Vel5','Pressure5','Bed5','Surface5','Thickness5','SmbMassBalance5'};
field_tolerances={1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,...
	1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,...
	1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,...
	1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,...
	1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13};
field_values={...
	(md.results.TransientSolution(400).Vx),...
	(md.results.TransientSolution(400).Vy),...
	(md.results.TransientSolution(400).Vz),...
	(md.results.TransientSolution(400).Vel),...
	(md.results.TransientSolution(400).Pressure),...
	(md.results.TransientSolution(400).Base),...
	(md.results.TransientSolution(400).Surface),...
	(md.results.TransientSolution(400).Thickness),...
	(md.results.TransientSolution(400).SmbMassBalance),...
	(md.results.TransientSolution(800).Vx),...
	(md.results.TransientSolution(800).Vy),...
	(md.results.TransientSolution(800).Vz),...
	(md.results.TransientSolution(800).Vel),...
	(md.results.TransientSolution(800).Pressure),...
	(md.results.TransientSolution(800).Base),...
	(md.results.TransientSolution(800).Surface),...
	(md.results.TransientSolution(800).Thickness),...
	(md.results.TransientSolution(800).SmbMassBalance),...
	(md.results.TransientSolution(1200).Vx),...
	(md.results.TransientSolution(1200).Vy),...
	(md.results.TransientSolution(1200).Vz),...
	(md.results.TransientSolution(1200).Vel),...
	(md.results.TransientSolution(1200).Pressure),...
	(md.results.TransientSolution(1200).Base),...
	(md.results.TransientSolution(1200).Surface),...
	(md.results.TransientSolution(1200).Thickness),...
	(md.results.TransientSolution(1200).SmbMassBalance),...
	(md.results.TransientSolution(1600).Vx),...
	(md.results.TransientSolution(1600).Vy),...
	(md.results.TransientSolution(1600).Vz),...
	(md.results.TransientSolution(1600).Vel),...
	(md.results.TransientSolution(1600).Pressure),...
	(md.results.TransientSolution(1600).Base),...
	(md.results.TransientSolution(1600).Surface),...
	(md.results.TransientSolution(1600).Thickness),...
	(md.results.TransientSolution(1600).SmbMassBalance),...
	(md.results.TransientSolution(2000).Vx),...
	(md.results.TransientSolution(2000).Vy),...
	(md.results.TransientSolution(2000).Vz),...
	(md.results.TransientSolution(2000).Vel),...
	(md.results.TransientSolution(2000).Pressure),...
	(md.results.TransientSolution(2000).Base),...
	(md.results.TransientSolution(2000).Surface),...
	(md.results.TransientSolution(2000).Thickness),...
	(md.results.TransientSolution(2000).SmbMassBalance),...
	};

if printingflag,

	starttime = 360;
	endtime = 2000;
	res = 40;
	ts = [starttime:res:endtime];

	index = md.mesh.elements;
	x1=md.mesh.x(index(:,1)); x2=md.mesh.x(index(:,2)); x3=md.mesh.x(index(:,3));
	y1=md.mesh.y(index(:,1)); y2=md.mesh.y(index(:,2)); y3=md.mesh.y(index(:,3));
	areas=(0.5*((x2-x1).*(y3-y1)-(y2-y1).*(x3-x1)));

	thickness = [];
	volume = [];
	massbal = [];
	velocity = [];
	for t=starttime:endtime
		thickness = [thickness (md.results.TransientSolution(t).Thickness)];
		volume = [volume mean(md.results.TransientSolution(t).Thickness.value,2).*areas];
		massbal = [massbal (md.results.TransientSolution(t).SmbMassBalance)];
		velocity = [velocity (md.results.TransientSolution(t).Vel)];
	end

	figure('Position', [0 0 1060 1060])

	options = plotoptions('data','transient_movie','unit','km');
	options = options.list{1};
	options = checkplotoptions(md,options);

	%loop over the time steps
	results=md.results.TransientSolution;
	count = 1;
	for i=ts

		subplot(5,9,[28:31 37:40])
		set(gca,'pos',get(gca,'pos')+[-0.08 -0.08 0.07 0.08])
		field = 'Thickness';

		%process data
		[x y z elements is2d isplanet]=processmesh(md,results(i).(field),options);
		[data datatype]=processdata(md,results(i).(field),options);

		titlestring=[field ' at time ' num2str(results(i).time/md.constants.yts) ' year'];
		plot_unit(x,y,z,elements,data,is2d,isplanet,datatype,options)
		options=changefieldvalue(options,'title',titlestring);
		options=addfielddefault(options,'colorbar',1);
		options=changefieldvalue(options,'caxis',[0 max(max(thickness))]);
		applyoptions(md,[],options);

		subplot(5,9,[33:36 42:45])
		set(gca,'pos',get(gca,'pos')+[-0.01 -0.08 0.07 0.08])
		field = 'Vel';

		%process data
		[x y z elements is2d isplanet]=processmesh(md,results(i).(field),options);
		[data datatype]=processdata(md,results(i).(field),options);

		titlestring=[field ' at time ' num2str(results(i).time/md.constants.yts) ' year'];
		plot_unit(x,y,z,elements,data,is2d,isplanet,datatype,options)
		options=changefieldvalue(options,'title',titlestring);
		options=addfielddefault(options,'colorbar',1);
		options=changefieldvalue(options,'caxis',[0 max(max(velocity))]);
		applyoptions(md,[],options);

		subplot(5,4,1:4)
		cla
		set(gca,'pos',get(gca,'pos')+[-0.07 0.03 0.12 0.015])
		plot(starttime:endtime,mean(massbal),'k','LineWidth', 4)
		hold on
		ya = ylim;
		plot([i i], ya, 'r', 'LineWidth',6)
		ylim(ya); xlim([starttime endtime]);
		title('Surface Mass Balance','FontSize',14)
		ylabel('m/year','FontSize',14)

		subplot(5,4,5:8)
		cla
		set(gca,'pos',get(gca,'pos')+[-0.07 0.015 0.12 0.015])
		plot(starttime:endtime,sum(volume)/1000/1000/1000,'LineWidth',4)
		hold on
		ya = ylim;
		plot([i i], ya, 'r', 'LineWidth',6)
		ylim(ya); xlim([starttime endtime]);
		title('Ice Volume','FontSize',14)
		ylabel('km^3','FontSize',14)

		subplot(5,4,9:12)
		cla
		set(gca,'pos',get(gca,'pos')+[-0.07 0 0.12 0.015])
		plot(starttime:endtime,mean(velocity)/1000, 'LineWidth', 4)
		hold on
		ya = ylim;
		plot([i i], ya, 'r', 'LineWidth',6)
		ylim(ya); xlim([starttime endtime]);
		title('Mean Velocity','FontSize', 14)
		ylabel('km/year','FontSize', 14)
		xlabel('year','FontSize', 14)

		set(gcf,'Renderer','zbuffer','color','white'); %fixes a bug on Mac OS X (not needed in future Matlab version)
		if i==starttime,
			%initialize images and frame
			frame=getframe(gcf);
			[images,map]=rgb2ind(frame.cdata,256,'nodither');
			images(1,1,1,length(ts))=0;
		else
			frame=getframe(gcf);
			images(:,:,1,count) = rgb2ind(frame.cdata,map,'nodither');
		end

		count = count+1;

	end

	filename='transawtooth3d.gif';
	imwrite(images,map,filename,'DelayTime',1.0,'LoopCount',inf)

end %printingflag
