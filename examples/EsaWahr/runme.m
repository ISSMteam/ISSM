clear all;
addpath('../Functions');

steps=[1];

if any(steps==1)
	disp('   Step 1: Mesh creation');

	md=roundmesh(model,100000,10000);  % Domain radius and element size [m]

	save ./Models/EsaWahr_Mesh md;

	plotmodel(md,'data','mesh');
end

if any(steps==2)
	disp('   Step 2: Anisotropic mesh creation');

	md=roundmesh(model,100000,1000);

	disc_radius=20*1000;
	rad_dist=sqrt(md.mesh.x.^2+md.mesh.y.^2);	
	field = abs(rad_dist-disc_radius);

	md = bamg(md,'field',field,'err',50,'hmax',10000);

	save ./Models/EsaWahr_Mesh md;

	plotmodel (md,'data','mesh');
end

if any(steps==3)
	disp('   Step 3: Define loads');
	md = loadmodel('./Models/EsaWahr_Mesh');

	rho_w_i=md.materials.rho_freshwater/md.materials.rho_ice;

	index=md.mesh.elements;		
	x_cent=mean(md.mesh.x(index),2);
	y_cent=mean(md.mesh.y(index),2);

	md.esa.deltathickness = zeros(md.mesh.numberofelements,1);
	disc_radius=20; % [km]
	rad_dist=sqrt(x_cent.^2+y_cent.^2)/1000;	
	md.esa.deltathickness(rad_dist<=disc_radius) = -1.0*rho_w_i;

	save ./Models/EsaWahr_Loads md;

	plotmodel (md,'data',md.esa.deltathickness,'title','Ice height equivalent [m]');
end

if any(steps==4)
	disp('   Step 4: Parameterization');
	md = loadmodel('./Models/EsaWahr_Loads');

	love_numbers = lovenumbers('maxdeg',10000,'referenceframe','CF');
	md.esa.love_h = love_numbers.h;
	md.esa.love_l = love_numbers.l;

	md.mask.ice_levelset = -ones(md.mesh.numberofvertices,1);
	md.mask.ocean_levelset = ones(md.mesh.numberofvertices,1);

	di=md.materials.rho_ice/md.materials.rho_water;
	md.geometry.thickness=ones(md.mesh.numberofvertices,1);
	md.geometry.surface=(1-di)*zeros(md.mesh.numberofvertices,1);
	md.geometry.base=md.geometry.surface-md.geometry.thickness;
	md.geometry.bed=md.geometry.base;

	md.initialization.temperature=273.25*ones(md.mesh.numberofvertices,1);
	md.materials.rheology_B=paterson(md.initialization.temperature);
	md.materials.rheology_n=3*ones(md.mesh.numberofelements,1);

	md.miscellaneous.name='EsaWahr';

	save ./Models/EsaWahr_Parameterization md;
end

if any(steps==5)
	disp('   Step 5: Solve Esa solver');
	md = loadmodel('./Models/EsaWahr_Parameterization');

	md.esa.requested_outputs = {'EsaUmotion','EsaXmotion','EsaYmotion'};

	md.cluster=generic('name',oshostname(),'np',3);
	md.verbose=verbose('111111111');

	md=solve(md,'Esa');

	save ./Models/EsaWahr_Solution md;
end

if any(steps==6)
	disp('   Step 6: Plot solutions');
	md = loadmodel('./Models/EsaWahr_Solution');

	vert = md.results.EsaSolution.EsaUmotion*1000;		% [mm]
	horz_n = md.results.EsaSolution.EsaYmotion*1000;	% [mm]
	horz_e = md.results.EsaSolution.EsaXmotion*1000;	% [mm]
	horz = sqrt(horz_n.^2+horz_e.^2);						% [mm]

	set(0,'DefaultAxesFontSize',24,'DefaultAxesLineWidth',1,'DefaultTextFontSize',24,'DefaultLineMarkerSize',6);
	figure('Position', [100, 100, 800, 600]);
	plotmodel(md,'data',vert,...
		'xTickLabel#all',[],'yTickLabel#all',[],...
		'colormap#all','jet','colorbarcornerposition#all','south',...
		'expdisp#all','./RoundDomain.exp',...
		'gridded#all',1,...
		'axispos#1',[0.02 0.505 0.47 0.47],...
		'colorbarpos#1',[0.045,0.55,0.18,0.02],'colorbartitle#1','Vertical [mm]',...
		'caxis#1',[0 3.5],...
		'data',horz,...
		'axispos#2',[0.505 0.505 0.47 0.47],...
		'colorbarpos',[0.53,0.55,0.18,0.02],'colorbartitle#2','Horizontal [mm]',...
		'caxis#2',[0 0.5],...
		'data',horz_n,...
		'axispos',[0.02 0.02 0.47 0.47],...
		'colorbarpos',[0.045,0.065,0.18,0.02],'colorbartitle#3','North-south [mm]',...
		'data',horz_e,...
		'caxis#3-4',[-0.5 0.5],...
		'axispos',[0.505 0.02 0.47 0.47],...
		'colorbarpos',[0.53,0.065,0.18,0.02],'colorbartitle#4','East-west [mm]');
	%export_fig('Fig5.pdf');
end

if any(steps==7)
	disp('   Step 7: Compare results against Wahr semi-analytic solutions');
	md = loadmodel('./Models/EsaWahr_Solution');

	vert = md.results.EsaSolution.EsaUmotion*1000;		% [mm]
	horz_n = md.results.EsaSolution.EsaYmotion*1000;	% [mm]
	horz_e = md.results.EsaSolution.EsaXmotion*1000;	% [mm]
	horz = sqrt(horz_n.^2+horz_e.^2);						% [mm]

	xi=[0:500:100000]; % grid points [m]
	yi=zeros(1,length(xi));
	vert_track=griddata(md.mesh.x,md.mesh.y,vert,xi,yi,'linear');
	horz_track=griddata(md.mesh.x,md.mesh.y,horz,xi,yi,'linear');

	% semi-analytic solution (Wahr et al., 2013, JGR, Figure 1)
	disc_radius = 20*1000; % [m]
	[vert_wahr, horz_wahr, acc_wahr] = wahr(disc_radius,xi,md.esa.love_h,md.esa.love_l,md.esa.love_h);
	vert_wahr = vert_wahr*1e3; % mm 
	horz_wahr = horz_wahr*1e3; % mm 

	set(0,'DefaultAxesFontSize',16,'DefaultAxesLineWidth',1,'DefaultTextFontSize',16,'DefaultLineMarkerSize',6);
	figure1=figure('Position', [100, 100, 700, 400]);
	ylabel_1={'0',' ','1','','2','','3',''};
	axes1 = axes('Layer','top','Position',[0.1 0.15 0.8 0.8],...
		'XTick',[0:10:100],'xlim',[0 100],...
		'ylim',[0 3.5],'Ytick',[0:0.5:3.5],'yticklabel',ylabel_1);
		box(axes1,'on'); hold(axes1,'all'); grid on;
		xlabel(axes1,'Radial distance [km]');
		ylabel(axes1,'Displacement [mm]');
		plot([20 20],[0 3.5],'-k','linewidth',2,'parent',axes1);
		% analytic soln
		h3=plot(xi/1000,vert_wahr,'-r','linewidth',5,'parent',axes1);
		h4=plot(xi/1000,horz_wahr,'-m','linewidth',5,'parent',axes1);
		% ISSM soln
		h1=plot(xi/1000,vert_track*917/1000,'-b','linewidth',3,'parent',axes1);
		h2=plot(xi/1000,horz_track*917/1000,'-c','linewidth',3,'parent',axes1);
		ag1 = gca;
		leg1a = legend(ag1,[h3,h1,h4,h2],'Vertical (Wahr)','Vertical (ISSM)','Horizontal (Wahr)','Horizontal (ISSM)');
		set(leg1a,'location','east','Orientation','Vertical','Box','Off','FontSize',16);
		set(gcf,'color','w');
	%export_fig('Fig6.pdf');
end
