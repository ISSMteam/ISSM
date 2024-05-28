%clear all;
addpath('../Data','../Functions');

%steps=[1:7];

if any(steps==1) % {{{ 
	disp('   Step 1: Global mesh creation');

	numrefine=1;
	resolution=150*1e3;			% inital resolution [m]
	radius = 6.371012*10^6;		% mean radius of Earth, m
	mindistance_coast=150*1e3;	% coastal resolution [m]
	mindistance_land=300*1e3;	% resolution on the continents [m]
	maxdistance=600*1e3;		% max element size (on mid-oceans) [m]

	md=model;
	md.mesh=gmshplanet('radius',radius*1e-3,'resolution',resolution*1e-3); % attributes in [km]

	for i=1:numrefine,

		ocean_mask=gmtmask(md.mesh.lat,md.mesh.long);

		distance=zeros(md.mesh.numberofvertices,1);

		pos=find(~ocean_mask);	coaste.lat=md.mesh.lat(pos);	coaste.long=md.mesh.long(pos);		pos=find(ocean_mask);	coasto.lat=md.mesh.lat(pos);	coasto.long=md.mesh.long(pos);
		for j=1:md.mesh.numberofvertices
			phi1=md.mesh.lat(j)/180*pi; lambda1=md.mesh.long(j)/180*pi;
			if ocean_mask(j),
				phi2=coaste.lat/180*pi; lambda2=coaste.long/180*pi;
				deltaphi=abs(phi2-phi1); deltalambda=abs(lambda2-lambda1);
				d=radius*2*asin(sqrt(sin(deltaphi/2).^2+cos(phi1).*cos(phi2).*sin(deltalambda/2).^2));
			else
				phi2=coasto.lat/180*pi; lambda2=coasto.long/180*pi;
				deltaphi=abs(phi2-phi1); deltalambda=abs(lambda2-lambda1);
				d=radius*2*asin(sqrt(sin(deltaphi/2).^2+cos(phi1).*cos(phi2).*sin(deltalambda/2).^2));
			end
			distance(j)=min(d);
		end
		pos=find(distance<mindistance_coast); distance(pos)=mindistance_coast;

		pos2=find(ocean_mask~=1 & distance>mindistance_land);
		distance(pos2)=mindistance_land;

		dist=min(maxdistance,distance);
		md.mesh=gmshplanet('radius',radius*1e-3,'resolution',resolution*1e-3,'refine',md.mesh,'refinemetric',dist);
	end

	ocean_mask=gmtmask(md.mesh.lat,md.mesh.long);
	pos = find(ocean_mask==0);	md.mask.ocean_levelset=-ones(md.mesh.numberofvertices,1);	md.mask.ocean_levelset(pos)=1;
	save ./Models/SlrGRACE_Mesh md;

	plotmodel (md,'data',md.mask.ocean_levelset,'edgecolor','k','view',[45 45]);

end % }}} 
if any(steps==2) % {{{ 
	disp('   Step 2: Define loads in meters of ice height equivalent');
	md = loadmodel('./Models/SlrGRACE_Mesh');

	year_month = 2007+15/365;
	time_range = [year_month year_month];

	onvertex = 1; % map data on vertex. If 0, it maps on the elemental centroid. 
	water_load = grace(md.mesh.elements,md.mesh.lat,md.mesh.long,time_range(1),time_range(2),onvertex);
	rho_water2ice = md.materials.rho_freshwater/md.materials.rho_ice; 
	ice_load = water_load*rho_water2ice; % ice height equivalent. 

	%Geometry for the bed, arbitrary thickness of 100:
	md.geometry.bed=zeros(md.mesh.numberofvertices,1);
	md.geometry.base=md.geometry.bed;
	md.geometry.thickness = 100*ones(md.mesh.numberofvertices,1);
	md.geometry.surface=md.geometry.bed+md.geometry.thickness;

	md.masstransport.spcthickness = repmat(md.geometry.thickness,1,2); 
	md.masstransport.spcthickness(:,end) = md.masstransport.spcthickness(:,end) + ice_load; 
	md.masstransport.spcthickness(end+1,:) = [0 1]; % dummy time. 

	md.smb.mass_balance=zeros(md.mesh.numberofvertices,1); 

	save ./Models/SlrGRACE_Loads md;

	plotmodel (md,'data',ice_load,...
		'edgecolor','k','view',[45 45],'caxis',[-.1 .1],...
		'title','Ice height equivalent [m]');
end % }}} 
if any(steps==3) % {{{ 
	disp('   Step 3: Parameterization');
	md = loadmodel('./Models/SlrGRACE_Loads');
	
	md.mask.ice_levelset=-md.mask.ocean_levelset;
	
	md.solidearth.lovenumbers = lovenumbers('maxdeg',10000);
	md.solidearth.settings.reltol=NaN;
	md.solidearth.settings.abstol=1e-3;
	md.solidearth.settings.sealevelloading=1;
	md.solidearth.settings.isgrd=1;
	md.solidearth.settings.grdmodel=1;
	md.solidearth.settings.maxiter=10; 
	
	%time stepping:
	md.timestepping.start_time=0;
	md.timestepping.time_step=1;
	md.timestepping.final_time=1;

	%masstransport:
	md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
	md.basalforcings.floatingice_melting_rate=zeros(md.mesh.numberofvertices,1);
	md.initialization.vx=zeros(md.mesh.numberofvertices,1);
	md.initialization.vy=zeros(md.mesh.numberofvertices,1);
	md.initialization.sealevel=zeros(md.mesh.numberofvertices,1);
	md.initialization.str=0;

	md.miscellaneous.name='SlrGRACE';

	save ./Models/SlrGRACE_Parameterization md;
end % }}} 
if any(steps==4) % {{{ 
	disp('   Step 4: Solve Slr solver');
	md = loadmodel('./Models/SlrGRACE_Parameterization');

	md.solidearth.settings.viscous=0;
	md.solidearth.settings.selfattraction=1;
	md.solidearth.settings.elastic=1;
	md.solidearth.settings.rotation=1;

	%Physics:
	md.transient.issmb=0;
	md.transient.isstressbalance=0;
	md.transient.isthermal=0;
	md.transient.ismasstransport=1;
	md.transient.isslc=1;
	
	md.solidearth.requested_outputs={'Sealevel','Bed'}; 

	md.cluster=generic('name',oshostname(),'np',3);
	md.verbose=verbose('111111111');

	md=solve(md,'Transient');

	save ./Models/SlrGRACE_Solution md;
end % }}}
if any(steps==5) % {{{ 
	disp('   Step 5: Plot solutions');
	md = loadmodel('./Models/SlrGRACE_Solution');

	sol1 = (md.masstransport.spcthickness(1:end-1,2)-md.masstransport.spcthickness(1:end-1,1))*100; % [cm]
	sol2 = (md.results.TransientSolution.Sealevel-md.results.TransientSolution.Bed)*1000;	% [mm]

	sol_name={'Change in water equivalent height [cm]', 'Relative sea level [mm]'};
	fig_name={'Fig_dH.pdf','Fig_slr.pdf'};

	res = 1.0;

	[lat_grid, lon_grid] = meshgrid(linspace(-90,90,180/res), linspace(-180,180,360/res));
	sol_grid = zeros(size(lat_grid));

	for kk=1:2
		sol=eval(sprintf('sol%d',kk));

		if length(sol)==md.mesh.numberofelements
			for jj=1:md.mesh.numberofelements
				ii=(jj-1)*3;
				pp(ii+1:ii+3)=md.mesh.elements(jj,:);
			end
			for jj=1:md.mesh.numberofvertices
				pos=ceil(find(pp==jj)/3);
				temp(jj)=mean(sol(pos));
			end
			sol=temp';
		end

		F = scatteredInterpolant(md.mesh.lat,md.mesh.long,sol);
		F.Method = 'linear';
		F.ExtrapolationMethod = 'linear';

		sol_grid = F(lat_grid, lon_grid);
		sol_grid(isnan(sol_grid))=0;
		sol_grid(lat_grid>85 & sol_grid==0) = NaN;

		set(0,'DefaultAxesFontSize',18,'DefaultAxesLineWidth',1,'DefaultTextFontSize',18,'DefaultLineMarkerSize',8)
		figure1=figure('Position', [100, 100, 1000, 500]);
		gcf; load coastlines; cla;
		pcolor(lon_grid,lat_grid,sol_grid); shading flat; hold on;
		if (kk==1)
			geoshow(flipud(coastlat),flipud(coastlon),'DisplayType','polygon','FaceColor','white');
		else
			geoshow(coastlat,coastlon,'DisplayType','polygon','FaceColor','none');
		end
		plot(coastlon, coastlat,'k'); hold off;
		c1=colorbar;
		colormap('haxby');
		caxis([-floor(min(abs(min(sol)),abs(max(sol)))) floor(min(abs(min(sol)),abs(max(sol))))]);
		xlim([-180 180]);
		ylim([-90 90]);
		grid on;
		title(sol_name(kk));
		set(gcf,'color','w');
		%export_fig(fig_name{kk});
	end
end % }}} 
if any(steps==6) % {{{
	disp('   Step 6: Transient run');
	md = loadmodel('./Models/SlrGRACE_Parameterization');

	disp('Projecting  loads onto the mesh...');
	time_range = 2007 + [15 350]/365;
	onvertex = 1; % map data on vertex. If 0, it maps on the elemental centroid. 
	water_load = grace(md.mesh.elements,md.mesh.lat,md.mesh.long,time_range(1),time_range(2),onvertex);
	rho_water2ice = md.materials.rho_freshwater/md.materials.rho_ice; 
	ice_load = water_load*rho_water2ice; % ice height equivalent. 

	% masstransport evalulates diff between the successive times, so we should cumsum. 
	num_time = size(ice_load,2); 
	md.masstransport.spcthickness = repmat(md.geometry.thickness,1,num_time+1); 
	md.masstransport.spcthickness(:,2:end) = md.masstransport.spcthickness(:,2:end) + ice_load;
	md.masstransport.spcthickness(end+1,:) = 0:num_time; % dummy time. 

	%Physics 
	md.transient.issmb=0;
	md.transient.isstressbalance=0;
	md.transient.isthermal=0;
	md.transient.ismasstransport=1;
	md.transient.isslc=1;

	md.solidearth.settings.viscous=0;
	md.solidearth.settings.selfattraction=1;
	md.solidearth.settings.elastic=1;
	md.solidearth.settings.rotation=1;
	
	%time stepping:
	md.timestepping.start_time=0;
	md.timestepping.time_step=1;
	md.timestepping.final_time=num_time; 

	md.cluster=generic('name',oshostname(),'np',3);
	md.verbose=verbose('111111111');

	md.solidearth.requested_outputs = {'Sealevel','Bed'};

	md=solve(md,'tr');

	save ./Models/SlrGRACE_Transient md;
end % }}} 
if any(steps==7) % {{{ 
	disp('   Step 7: Plot transient');
	md = loadmodel('./Models/SlrGRACE_Transient');

	time = md.masstransport.spcthickness(end,2:end);

	for tt=1:length(time)
		gmsl(tt) = md.results.TransientSolution(tt).CumBslc*1000; % [mm]
		sol1(:,tt) = (md.masstransport.spcthickness(1:end-1,tt+1)-md.geometry.thickness)*100; % [cm]
		sol2(:,tt) = (md.results.TransientSolution(tt).Sealevel-md.results.TransientSolution(tt).Bed)*1000;	% [mm]
	end
	sol_name = {'Change in water equivalent height [cm]', 'Relative sea level [mm]'};
	movie_name = {'Movie_dH','Movie_slr'};
	
	res = 1.0;

	[lat_grid, lon_grid] = meshgrid(linspace(-90,90,180/res), linspace(-180,180,360/res));
	sol_grid = zeros(size(lat_grid));

	for kk=1:2
		sol=eval(sprintf('sol%d',kk));

		if length(sol)==md.mesh.numberofelements
			for jj=1:md.mesh.numberofelements
				ii=(jj-1)*3;
				pp(ii+1:ii+3)=md.mesh.elements(jj,:);
			end
			for jj=1:md.mesh.numberofvertices
				pos=ceil(find(pp==jj)/3);
				temp2(jj,:)=mean(sol(pos,:));
			end
			sol=temp2;
		end

		vidObj = VideoWriter(movie_name{kk});
		vidObj.FrameRate=2; % frames per second
		open(vidObj);

		for jj=1:length(time)
			fprintf('creating frame %d...', jj);
			F = scatteredInterpolant(md.mesh.lat,md.mesh.long,sol(:,jj));
			F.Method = 'linear';
			F.ExtrapolationMethod = 'linear';

			sol_grid = F(lat_grid, lon_grid);
			sol_grid(isnan(sol_grid))=0;
			sol_grid(lat_grid>85 & sol_grid==0) = NaN;

			set(0,'DefaultAxesFontSize',18,'DefaultAxesLineWidth',1,'DefaultTextFontSize',18,'DefaultLineMarkerSize',8)
			figure1=figure('Position', [100, 100, 1000, 500]);
			gcf; load coastlines; cla;
			pcolor(lon_grid,lat_grid,sol_grid); shading flat; hold on;
			if (kk==1)
				geoshow(flipud(coastlat),flipud(coastlon),'DisplayType','polygon','FaceColor','white');
			else
				geoshow(coastlat,coastlon,'DisplayType','polygon','FaceColor','none');
			end
			plot(coastlon,coastlat,'k'); hold off;
			c1=colorbar;
			colormap('haxby');
			%caxis([-floor(min(abs(min(min(sol))),abs(max(max(sol))))) floor(min(abs(min(min(sol))),abs(max(max(sol)))))]);
			xlim([-180 180]);
			ylim([-90 90]);
			grid on;
			title(sol_name(kk));
			set(gcf,'color','w');
			writeVideo(vidObj,getframe(gcf));
			close
			fprintf('done\n');
		end
		disp('closing vidObj...');
		close(vidObj);
	end

	% plot GMSL time series
	plot(time,gmsl,'-*','linewidth',3); grid on;
	xlabel('# month');
	ylabel('GMSL [mm/yr]');
	set(gcf,'color','w');
	%export_fig('Fig7.pdf');
end % }}} 

