clear all;

steps=[4];

if any(steps==1) 
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
		ocean_mask_levelset=gmtmask(md.mesh.lat,md.mesh.long);

		distance=zeros(md.mesh.numberofvertices,1);

		pos=find(~ocean_mask_levelset);	coaste.lat=md.mesh.lat(pos);	coaste.long=md.mesh.long(pos);
		pos=find(ocean_mask_levelset);	coasto.lat=md.mesh.lat(pos);	coasto.long=md.mesh.long(pos);

		for j=1:md.mesh.numberofvertices
			phi1=md.mesh.lat(j)/180*pi; lambda1=md.mesh.long(j)/180*pi;
			if ocean_mask_levelset(j),
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

		pos2=find(ocean_mask_levelset~=1 & distance>mindistance_land);
		distance(pos2)=mindistance_land;

		dist=min(maxdistance,distance);
		md.mesh=gmshplanet('radius',radius*1e-3,'resolution',resolution*1e-3,'refine',md.mesh,'refinemetric',dist);
	end

	ocean_mask=gmtmask(md.mesh.lat,md.mesh.long);
	pos = find(ocean_mask==0); 
	md.mask.ocean_levelset=-ones(md.mesh.numberofvertices,1); 
	md.mask.ocean_levelset(pos)=1; 

	save ./Models/SlrFarrell_Mesh md;

	plotmodel (md,'data',md.mask.ocean_levelset,'edgecolor','k');
end 

if any(steps==2) 
	disp('   Step 2: Define source as in Farrell, 1972, Figure 1');
	md = loadmodel('./Models/SlrFarrell_Mesh');

	md.solidearth.surfaceload.icethicknesschange=zeros(md.mesh.numberofvertices,1);
	md.mask.ice_levelset=ones(md.mesh.numberofvertices,1); 

	pos = find(md.mask.ocean_levelset==-1); 
	md.solidearth.initialsealevel=zeros(md.mesh.numberofvertices,1);
	md.solidearth.initialsealevel(pos)=1; % 1 m SLR everywhere
	md.dsl.global_average_thermosteric_sea_level_change=[0;0];
	md.dsl.sea_surface_height_change_above_geoid=zeros(md.mesh.numberofvertices+1,1);
	md.dsl.sea_water_pressure_change_at_sea_floor=zeros(md.mesh.numberofvertices+1,1);

	save ./Models/SlrFarrell_Loads md;

	plotmodel (md,'data',md.solidearth.initialsealevel,'view',[90 90],'title#all','Initial sea-level [m]');
end 

if any(steps==3) 
	disp('   Step 3: Parameterization');
	md = loadmodel('./Models/SlrFarrell_Loads');

	md.solidearth.lovenumbers=lovenumbers('maxdeg',10000);

	%md.mask.ice_levelset = ones(md.mesh.numberofvertices,1);
	%md.mask.ocean_levelset = -ones(md.mesh.numberofvertices,1);
	%pos=find(md.mesh.lat <-80);
	%md.mask.ice_levelset(pos(1))=-1; % ice yes!
	%md.mask.ocean_levelset(pos(1))=1; % ice grounded!

	% arbitary to pass consistency check. 
	md.geometry.bed=-ones(md.mesh.numberofvertices,1);
	md.geometry.surface=ones(md.mesh.numberofvertices,1);
	md.geometry.base=md.geometry.bed; 
	md.geometry.thickness=md.geometry.surface-md.geometry.base; 

	md.initialization.temperature=273.25*ones(md.mesh.numberofvertices,1);
	md.materials.rheology_B=paterson(md.initialization.temperature);
	md.materials.rheology_n=3*ones(md.mesh.numberofelements,1);

	md.miscellaneous.name='SlrFarrell';

	save ./Models/SlrFarrell_Parameterization md;
end 

if any(steps==4) 
	disp('   Step 4: Solve Slr solver');
	md = loadmodel('./Models/SlrFarrell_Parameterization');

	md.cluster=generic('name',oshostname(),'np',3);
	md.verbose=verbose('111111111');

	md.solidearth.settings.reltol = 0.1/100; % percent change in solution

	md=solve(md,'Slr');

	save ./Models/SlrFarrell_Solution md;
end 

if any(steps==5) 
	disp('   Step 5: Plot solutions');
	md = loadmodel('./Models/SlrFarrell_Solution');

	sol = md.results.SealevelriseSolution.Sealevel*100; % per cent normalized by GMSL (which 1 m)

	res = 1; % [degree]

	[lat_grid, lon_grid] = meshgrid(linspace(-90,90,180/res), linspace(-180,180,360/res));
	sol_grid = zeros(size(lat_grid));

	F = scatteredInterpolant(md.mesh.lat,md.mesh.long,sol);
	F.Method = 'linear';
	F.ExtrapolationMethod = 'linear';

	sol_grid = F(lat_grid, lon_grid);
	sol_grid(isnan(sol_grid))=0;
	sol_grid(lat_grid>85 & sol_grid==0)=NaN;

	set(0,'DefaultAxesFontSize',18,'DefaultAxesLineWidth',1,'DefaultTextFontSize',18,'DefaultLineMarkerSize',8)
	figure1=figure('Position', [100, 100, 1000, 500]);
	gcf; load coast; cla;
	pcolor(lon_grid,lat_grid,sol_grid); shading flat; hold on;
	[C,h]=contour(lon_grid,lat_grid,sol_grid,[96 98 100 102 104 105],'-k','linewidth',2);
	clabel(C,h,'FontSize',18,'Color','red','LabelSpacing',500);
	geoshow(lat,long,'DisplayType','polygon','FaceColor',[.82 .82 .82]);
	plot(long,lat,'k'); hold off;
	c1=colorbar;
	colormap(flipud(haxby));
	caxis([96 105]);
	xlim([-170 170]);
	ylim([-85 85]);
	grid on;
	title('Relative sea-level [% of GMSL]');
	set(gcf,'color','w');
	%export_fig('Fig5.pdf');
end 
