clear all;
addpath('../Data','../Functions');

steps=[1];

if any(steps==1)
	disp('   Step 1: Global mesh creation');

	resolution=300;			% [km]
	radius = 6.371012*10^3;	% [km]

	md=model;
	md.mesh=gmshplanet('radius',radius,'resolution',resolution);

	md.mask.ocean_levelset=gmtmask(md.mesh.lat,md.mesh.long);

	save ./Models/EsaGRACE_Mesh md;

	plotmodel (md,'data',md.mask.ocean_levelset,'edgecolor','k');
end

if any(steps==2)
	disp('   Step 2: Define loads in meters of ice height equivalent');
	md = loadmodel('./Models/EsaGRACE_Mesh');

	year_month = 2007+15/365;
	time_range = [year_month year_month];

	onvertex = 0; % map data on vertex. If 0, it maps on the elemental centroid. 
	water_load = grace(md.mesh.elements,md.mesh.lat,md.mesh.long,time_range(1),time_range(2),onvertex);

	md.esa.deltathickness = water_load*md.materials.rho_freshwater/md.materials.rho_ice; % ice height equivalent

	save ./Models/EsaGRACE_Loads md;

	plotmodel (md,'data',md.esa.deltathickness,...
		'view',[90 -90],'caxis',[-.1 .1],...
		'title','Ice height equivalent [m]');
end

if any(steps==3)
	disp('   Step 3: Parameterization');
	md = loadmodel('./Models/EsaGRACE_Loads');

	love_numbers = lovenumbers('maxdeg',10000,'referenceframe','CF');
	md.esa.love_h = love_numbers.h;
	md.esa.love_l = love_numbers.l;

	md.mask.ocean_levelset = ones(md.mesh.numberofvertices,1);
	md.mask.ice_levelset = ones(md.mesh.numberofvertices,1);
	pos=find(md.esa.deltathickness~=0);
	md.mask.ice_levelset(md.mesh.elements(pos,:))=-1;
	%md.mask.land_levelset = 1-md.mask.ocean_levelset;

	di=md.materials.rho_ice/md.materials.rho_water;
	md.geometry.thickness=ones(md.mesh.numberofvertices,1);
	md.geometry.surface=(1-di)*zeros(md.mesh.numberofvertices,1);
	md.geometry.base=md.geometry.surface-md.geometry.thickness;
	md.geometry.bed=md.geometry.base;

	md.initialization.temperature=273.25*ones(md.mesh.numberofvertices,1);
	md.materials.rheology_B=paterson(md.initialization.temperature);
	md.materials.rheology_n=3*ones(md.mesh.numberofelements,1);

	md.miscellaneous.name='EsaGRACE';

	save ./Models/EsaGRACE_Parameterization md;
end

if any(steps==4)
	disp('   Step 4: Solve Esa solver');
	md = loadmodel('./Models/EsaGRACE_Parameterization');

	md.esa.requested_outputs = {'EsaUmotion','EsaNmotion','EsaEmotion'};

	md.cluster=generic('name',oshostname(),'np',3);
	md.verbose=verbose('111111111');

	md=solve(md,'Esa');

	save ./Models/EsaGRACE_Solution md;
end

if any(steps==5)
	disp('   Step 5: Plot solutions');
	md = loadmodel('./Models/EsaGRACE_Solution');

	sol1 = md.esa.deltathickness*100;				% [cm]
	sol2 = md.results.EsaSolution.EsaUmotion*1000;	% [mm]
	sol3 = md.results.EsaSolution.EsaNmotion*1000;	% [mm]
	sol4 = md.results.EsaSolution.EsaEmotion*1000;	% [mm]

	sol_name={'Change in water equivalent height [cm]', 'Vertical displacement [mm]',...
		'Horizontal (NS) displacement [mm]', 'Horizontal (EW) displacement [mm]'};
	fig_name={'Fig_dH.pdf','Fig_vert.pdf','Fig_horzNS.pdf','Fig_horzEW.pdf'};

	res = 1.0; % [degree]

	[lat_grid, lon_grid] = meshgrid(linspace(-90,90,180/res), linspace(-180,180,360/res));
	sol_grid = zeros(size(lat_grid));

	for kk=1:4
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
		sol_grid(lat_grid>85 & sol_grid==0)=NaN;

		set(0,'DefaultAxesFontSize',18,'DefaultAxesLineWidth',1,'DefaultTextFontSize',18,'DefaultLineMarkerSize',8)
		figure1=figure('Position', [100, 100, 1000, 500]);
		gcf; load coastlines; cla;
		pcolor(lon_grid,lat_grid,sol_grid); shading flat; hold on;
		if (kk==1)
			geoshow(flipud(coastlat),flipud(coastlon),'DisplayType','polygon','FaceColor','white');
		end
		plot(coastlon,coastlat,'k'); hold off;
		c1=colorbar;
		colormap('haxby');
		caxis([-min(abs(min(sol)),abs(max(sol))) min(abs(min(sol)),abs(max(sol)))]);
		xlim([-180 180]);
		ylim([-90 90]);
		grid on;
		title(sol_name(kk));
		set(gcf,'color','w');
		%export_fig(fig_name{kk});
	end
end
