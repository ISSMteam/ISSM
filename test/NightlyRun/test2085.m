%Test Name: FourierLoveKernels
% Homogeneous Earth, for which analytic solutions exist. 
% Love kernels for degree 2, 20, 200 (tested against analytic solns).  
% skip benchmarking for the inner-most interface. 

% set validation=1 for comparing against the analytic solutions. 
validation=0; 

% for volumetric potential
md=model();
md.groundingline.migration='None';

md.materials=materials('litho');
cst=365.25*24*3600*1000;

md.materials.numlayers = 40;
md.love.forcing_type = 9;

md.materials.density=zeros(md.materials.numlayers,1)+5511;
md.materials.lame_mu=zeros(md.materials.numlayers,1)+0.75e11;
md.materials.lame_lambda=zeros(md.materials.numlayers,1)+5e17;
md.materials.issolid=ones(md.materials.numlayers,1);
md.materials.rheologymodel=zeros(md.materials.numlayers,1);

%the following isn't used here but needs to have arrays of consistent size with the rest of the materials
md.materials.viscosity=zeros(md.materials.numlayers,1)+1e21;
md.materials.burgers_mu=md.materials.lame_mu/3;
md.materials.burgers_viscosity=md.materials.viscosity/10;
md.materials.ebm_alpha= ones(md.materials.numlayers,1)*0.9;
md.materials.ebm_delta= ones(md.materials.numlayers,1)*0.2;
md.materials.ebm_taul= ones(md.materials.numlayers,1)*54*60; %54min
md.materials.ebm_tauh= ones(md.materials.numlayers,1)*18.6*cst/1e3; %18.6yr

md.materials.radius =  linspace(10e3,6371e3,md.materials.numlayers+1)';
md.love.g0 = 9.8134357285509388; % directly grabbed from fourierlovesolver for this particular case. 

md.love.allow_layer_deletion=1;
md.love.frequencies=0;
md.love.nfreq=1;
md.love.istemporal=0;

md.love.sh_nmin = 2;
md.love.sh_nmax = 200;
md.love.love_kernels=1; 

md.miscellaneous.name='kernels';
md.cluster=generic('name',oshostname(),'np',3);
md.verbose=verbose('111111101');

md=solve(md,'lv');

% save yi's for all layers except for the inner-most one, at select degrees. 
degrees = [2 20 200];  % we archive solutions for degrees 2, 20, 200 
kernels=reshape(md.results.LoveSolution.LoveKernels, [md.love.sh_nmax+1 md.love.nfreq md.materials.numlayers+1 6]);

% extract love kernels {{{ 
% degree 2.
y1_tidal_degree002 = squeeze(kernels(degrees(1)+1,1,2:end,1));
y2_tidal_degree002 = squeeze(kernels(degrees(1)+1,1,2:end,2));
y3_tidal_degree002 = squeeze(kernels(degrees(1)+1,1,2:end,3));
y4_tidal_degree002 = squeeze(kernels(degrees(1)+1,1,2:end,4));
y5_tidal_degree002 = squeeze(kernels(degrees(1)+1,1,2:end,5));
y6_tidal_degree002 = squeeze(kernels(degrees(1)+1,1,2:end,6));

% degree 20. 
y1_tidal_degree020 = squeeze(kernels(degrees(2)+1,1,2:end,1));
y2_tidal_degree020 = squeeze(kernels(degrees(2)+1,1,2:end,2));
y3_tidal_degree020 = squeeze(kernels(degrees(2)+1,1,2:end,3));
y4_tidal_degree020 = squeeze(kernels(degrees(2)+1,1,2:end,4));
y5_tidal_degree020 = squeeze(kernels(degrees(2)+1,1,2:end,5));
y6_tidal_degree020 = squeeze(kernels(degrees(2)+1,1,2:end,6));

% degree 200. 
y1_tidal_degree200 = squeeze(kernels(degrees(3)+1,1,2:end,1));
y2_tidal_degree200 = squeeze(kernels(degrees(3)+1,1,2:end,2));
y3_tidal_degree200 = squeeze(kernels(degrees(3)+1,1,2:end,3));
y4_tidal_degree200 = squeeze(kernels(degrees(3)+1,1,2:end,4));
y5_tidal_degree200 = squeeze(kernels(degrees(3)+1,1,2:end,5));
y6_tidal_degree200 = squeeze(kernels(degrees(3)+1,1,2:end,6));
% }}} 

% validate tidal potential solutions against the analytic solutions. {{{ 
if validation

	param.rho = md.materials.density(1); 
	param.mu = md.materials.lame_mu(1); 
	param.G = 6.67e-11;
	param.radius = md.materials.radius;  
	param.g0 = md.love.g0;
	param.source = md.love.forcing_type; 

	% check against analytic solutions at the following degrees 
	degrees_analytic = [2 4 8 16 32]; 
	% extract analytic solutions. 
	for jj=1:length(degrees_analytic)
		param.degree = degrees_analytic(jj); 
		y_temp = yi_analytic_homogenous(param); 
		if (exist('y_ana','var')~=1)
			y_ana = zeros([length(degrees_analytic) size(y_temp)]); 
		end
		y_ana(jj,:,:)=y_temp; 
	end
	y1_ana=squeeze(y_ana(:,:,1)); 
	y2_ana=squeeze(y_ana(:,:,2)); 
	y3_ana=squeeze(y_ana(:,:,3)); 
	y4_ana=squeeze(y_ana(:,:,4)); 
	y5_ana=squeeze(y_ana(:,:,5)); 
	y6_ana=squeeze(y_ana(:,:,6)); 

	depth = (max(param.radius)-param.radius)/1000; % km.

	kernels=reshape(md.results.LoveSolution.LoveKernels, [md.love.sh_nmax+1 md.love.nfreq md.materials.numlayers+1 6]);
	y1 = squeeze(kernels(:,1,:,1));
	y2 = squeeze(kernels(:,1,:,2));
	y3 = squeeze(kernels(:,1,:,3));
	y4 = squeeze(kernels(:,1,:,4));
	y5 = squeeze(kernels(:,1,:,5));
	y6 = squeeze(kernels(:,1,:,6));
	
	set(0,'DefaultAxesFontSize',16,'DefaultTextFontSize',15,'DefaultAxesLineWidth',1,...
		'DefaultTextInterpreter','Latex','DefaultAxesFontName','Arial','DefaultLineMarkerSize',8)
	figure1=figure('Position', [100, 100, 800, 1000]);
	%---------------------------------------------------------------------
	axes1 = axes('Layer','top','Position',[0.1 0.73 0.4 0.23]); 
	axis(axes1,[-0.004 0.088 -0.5 6.5]);
	ylabel(axes1,'Depth [$\times$1000 km]','fontsize',20); 
	xlabel(axes1,'$y_1$ [m]','fontsize',20); 
	box(axes1,'on'); grid(axes1,'on'); hold(axes1,'all');
	plot(y1_ana,depth/1000,'linewidth',2,'parent',axes1); 
	set(gca,'ColorOrderIndex',1,'Ydir','reverse'); 
	plot(y1(degrees_analytic+1,:),depth/1000,'o','linewidth',2,'parent',axes1); 
	%---------------------------------------------------------------------
	axes2 = axes('Layer','top','Position',[0.55 0.73 0.4 0.23]); 
	axis(axes2,[-200 3200 -0.5 6.5]);
	xlabel(axes2,'$y_2$ [Pa]','fontsize',20); 
	box(axes2,'on'); grid(axes2,'on'); hold(axes2,'all');
	plot(y2_ana,depth/1000,'linewidth',2,'parent',axes2); 
	set(gca,'ColorOrderIndex',1,'Ydir','reverse'); 
	plot(y2(degrees_analytic+1,:),depth/1000,'o','linewidth',2,'parent',axes2); 
	%---------------------------------------------------------------------
	axes3 = axes('Layer','top','Position',[0.1 0.4 0.4 0.23]); 
	axis(axes3,[-0.002 0.034 -0.5 6.5]);
	ylabel(axes3,'Depth [$\times$1000 km]','fontsize',20); 
	xlabel(axes3,'$y_3$ [m]','fontsize',20); 
	box(axes3,'on'); grid(axes3,'on'); hold(axes3,'all');
	plot(y3_ana,depth/1000,'linewidth',2,'parent',axes3); 
	set(gca,'ColorOrderIndex',1,'Ydir','reverse'); 
	plot(y3(degrees_analytic+1,:),depth/1000,'o','linewidth',2,'parent',axes3); 
	%---------------------------------------------------------------------
	axes4 = axes('Layer','top','Position',[0.55 0.4 0.4 0.23]); 
	axis(axes4,[-100 1600 -0.5 6.5]);
	xlabel(axes4,'$y_4$ [Pa]','fontsize',20); 
	box(axes4,'on'); grid(axes4,'on'); hold(axes4,'all');
	plot(y4_ana,depth/1000,'linewidth',2,'parent',axes4); 
	set(gca,'ColorOrderIndex',1,'Ydir','reverse'); 
	plot(y4(degrees_analytic+1,:),depth/1000,'o','linewidth',2,'parent',axes4); 
	%---------------------------------------------------------------------
	axes5 = axes('Layer','top','Position',[0.1 0.07 0.4 0.23]); 
	axis(axes5,[-0.05 1.5 0 6.5]);
	ylabel(axes5,'Depth [$\times$1000 km]','fontsize',20); 
	xlabel(axes5,'$y_5$ [m$^2$ s$^{-2}$]','fontsize',20); 
	box(axes5,'on'); grid(axes5,'on'); hold(axes5,'all');
	plot(y5_ana,depth/1000,'linewidth',2,'parent',axes5); 
	set(gca,'ColorOrderIndex',1,'Ydir','reverse'); 
	plot(y5(degrees_analytic+1,:),depth/1000,'o','linewidth',2,'parent',axes5); 
	%---------------------------------------------------------------------
	axes6 = axes('Layer','top','Position',[0.55 0.07 0.4 0.23]); 
	axis(axes6,[-1e-7 4e-7 0 6.5]);
	xlabel(axes6,'$y_6$ [m s$^{-2}$]','fontsize',20); 
	box(axes6,'on'); grid(axes6,'on'); hold(axes6,'all');
	plot(y6_ana,depth/1000,'linewidth',2,'parent',axes6); 
	set(gca,'ColorOrderIndex',1,'Ydir','reverse'); 
	plot(y6(degrees_analytic+1,:),depth/1000,'o','linewidth',2,'parent',axes6); 
	%---------------------------------------------------------------------

	legend(axes6,'n=2','n=4','n=8','n=16','n=32'); 
	%export_fig('/Users/adhikari/issm/trunk-jpl/src/m/contrib/adhikari/issm_vs_analytic_loading_homogeneous.pdf'); 
else
	% 
end 
% }}} 

% for surface load. 
md.love.forcing_type = 11;
md=solve(md,'lv');
kernels=reshape(md.results.LoveSolution.LoveKernels, [md.love.sh_nmax+1 md.love.nfreq md.materials.numlayers+1 6]);

% extract love kernels {{{ 
% degree 2. 
y1_loading_degree002 = squeeze(kernels(degrees(1)+1,1,2:end,1));
y2_loading_degree002 = squeeze(kernels(degrees(1)+1,1,2:end,2));
y3_loading_degree002 = squeeze(kernels(degrees(1)+1,1,2:end,3));
y4_loading_degree002 = squeeze(kernels(degrees(1)+1,1,2:end,4));
y5_loading_degree002 = squeeze(kernels(degrees(1)+1,1,2:end,5));
y6_loading_degree002 = squeeze(kernels(degrees(1)+1,1,2:end,6));

% degree 20. 
y1_loading_degree020 = squeeze(kernels(degrees(2)+1,1,2:end,1));
y2_loading_degree020 = squeeze(kernels(degrees(2)+1,1,2:end,2));
y3_loading_degree020 = squeeze(kernels(degrees(2)+1,1,2:end,3));
y4_loading_degree020 = squeeze(kernels(degrees(2)+1,1,2:end,4));
y5_loading_degree020 = squeeze(kernels(degrees(2)+1,1,2:end,5));
y6_loading_degree020 = squeeze(kernels(degrees(2)+1,1,2:end,6));

% degree 200. 
y1_loading_degree200 = squeeze(kernels(degrees(3)+1,1,2:end,1));
y2_loading_degree200 = squeeze(kernels(degrees(3)+1,1,2:end,2));
y3_loading_degree200 = squeeze(kernels(degrees(3)+1,1,2:end,3));
y4_loading_degree200 = squeeze(kernels(degrees(3)+1,1,2:end,4));
y5_loading_degree200 = squeeze(kernels(degrees(3)+1,1,2:end,5));
y6_loading_degree200 = squeeze(kernels(degrees(3)+1,1,2:end,6));
% }}} 

% validate loading solutions against the analytic solutions. {{{ 
if validation

	param.source = md.love.forcing_type; 

	% extract analytic solutions. 
	for jj=1:length(degrees_analytic)
		param.degree = degrees_analytic(jj); 
		y_temp = yi_analytic_homogenous(param); 
		if (exist('y_ana','var')~=1)
			y_ana = zeros([length(degrees_analytic) size(y_temp)]); 
		end
		y_ana(jj,:,:)=y_temp; 
	end
	y1_ana=squeeze(y_ana(:,:,1)); 
	y2_ana=squeeze(y_ana(:,:,2)); 
	y3_ana=squeeze(y_ana(:,:,3)); 
	y4_ana=squeeze(y_ana(:,:,4)); 
	y5_ana=squeeze(y_ana(:,:,5)); 
	y6_ana=squeeze(y_ana(:,:,6)); 

	depth = (max(param.radius)-param.radius)/1000; % km.

	kernels=reshape(md.results.LoveSolution.LoveKernels, [md.love.sh_nmax+1 md.love.nfreq md.materials.numlayers+1 6]);
	y1 = squeeze(kernels(:,1,:,1));
	y2 = squeeze(kernels(:,1,:,2));
	y3 = squeeze(kernels(:,1,:,3));
	y4 = squeeze(kernels(:,1,:,4));
	y5 = squeeze(kernels(:,1,:,5));
	y6 = squeeze(kernels(:,1,:,6));
	
	set(0,'DefaultAxesFontSize',16,'DefaultTextFontSize',15,'DefaultAxesLineWidth',1,...
		'DefaultTextInterpreter','Latex','DefaultAxesFontName','Arial','DefaultLineMarkerSize',8)
	figure1=figure('Position', [100, 100, 800, 1000]);
	%---------------------------------------------------------------------
	axes1 = axes('Layer','top','Position',[0.1 0.73 0.4 0.23]); 
	axis(axes1,[-0.11 0.01 -0.5 6.5]);
	ylabel(axes1,'Depth [$\times$1000 km]','fontsize',20); 
	xlabel(axes1,'$y_1$ [m]','fontsize',20); 
	box(axes1,'on'); grid(axes1,'on'); hold(axes1,'all');
	plot(y1_ana,depth/1000,'linewidth',2,'parent',axes1); 
	set(gca,'ColorOrderIndex',1,'Ydir','reverse'); 
	plot(y1(degrees_analytic+1,:),depth/1000,'o','linewidth',2,'parent',axes1); 
	%---------------------------------------------------------------------
	axes2 = axes('Layer','top','Position',[0.55 0.73 0.4 0.23]); 
	axis(axes2,[-1.8e4 0.2e4 -0.5 6.5]);
	xlabel(axes2,'$y_2$ [Pa]','fontsize',20); 
	box(axes2,'on'); grid(axes2,'on'); hold(axes2,'all');
	plot(y2_ana,depth/1000,'linewidth',2,'parent',axes2); 
	set(gca,'ColorOrderIndex',1,'Ydir','reverse'); 
	plot(y2(degrees_analytic+1,:),depth/1000,'o','linewidth',2,'parent',axes2); 
	%---------------------------------------------------------------------
	axes3 = axes('Layer','top','Position',[0.1 0.4 0.4 0.23]); 
	axis(axes3,[-0.023 0.002 -0.5 6.5]);
	ylabel(axes3,'Depth [$\times$1000 km]','fontsize',20); 
	xlabel(axes3,'$y_3$ [m]','fontsize',20); 
	box(axes3,'on'); grid(axes3,'on'); hold(axes3,'all');
	plot(y3_ana,depth/1000,'linewidth',2,'parent',axes3); 
	set(gca,'ColorOrderIndex',1,'Ydir','reverse'); 
	plot(y3(degrees_analytic+1,:),depth/1000,'o','linewidth',2,'parent',axes3); 
	%---------------------------------------------------------------------
	axes4 = axes('Layer','top','Position',[0.55 0.4 0.4 0.23]); 
	axis(axes4,[-1300 100 -0.5 6.5]);
	xlabel(axes4,'$y_4$ [Pa]','fontsize',20); 
	box(axes4,'on'); grid(axes4,'on'); hold(axes4,'all');
	plot(y4_ana,depth/1000,'linewidth',2,'parent',axes4); 
	set(gca,'ColorOrderIndex',1,'Ydir','reverse'); 
	plot(y4(degrees_analytic+1,:),depth/1000,'o','linewidth',2,'parent',axes4); 
	%---------------------------------------------------------------------
	axes5 = axes('Layer','top','Position',[0.1 0.07 0.4 0.23]); 
	axis(axes5,[-0.05 0.9 -0.5 6.5]);
	ylabel(axes5,'Depth [$\times$1000 km]','fontsize',20); 
	xlabel(axes5,'$y_5$ [m$^2$ s$^{-2}$]','fontsize',20); 
	box(axes5,'on'); grid(axes5,'on'); hold(axes5,'all');
	plot(y5_ana,depth/1000,'linewidth',2,'parent',axes5); 
	set(gca,'ColorOrderIndex',1,'Ydir','reverse'); 
	plot(y5(degrees_analytic+1,:),depth/1000,'o','linewidth',2,'parent',axes5); 
	%---------------------------------------------------------------------
	axes6 = axes('Layer','top','Position',[0.55 0.07 0.4 0.23]); 
	axis(axes6,[-1e-7 1e-6 -0.5 6.5]);
	xlabel(axes6,'$y_6$ [m s$^{-2}$]','fontsize',20); 
	box(axes6,'on'); grid(axes6,'on'); hold(axes6,'all');
	plot(y6_ana,depth/1000,'linewidth',2,'parent',axes6); 
	set(gca,'ColorOrderIndex',1,'Ydir','reverse'); 
	plot(y6(degrees_analytic+1,:),depth/1000,'o','linewidth',2,'parent',axes6); 
	%---------------------------------------------------------------------

	legend(axes6,'n=2','n=4','n=8','n=16','n=32'); 
	%export_fig('/Users/adhikari/issm/trunk-jpl/src/m/contrib/adhikari/issm_vs_analytic_tidal_homogeneous.pdf'); 
else
	% 
end 
% }}} 

field_names = {...
	'y1_tidal_degree002','y2_tidal_degree002','y3_tidal_degree002','y4_tidal_degree002','y5_tidal_degree002','y6_tidal_degree002',...
	'y1_tidal_degree020','y2_tidal_degree020','y3_tidal_degree020','y4_tidal_degree020','y5_tidal_degree020','y6_tidal_degree020',...
	'y1_tidal_degree200','y2_tidal_degree200','y3_tidal_degree200','y4_tidal_degree200','y5_tidal_degree200','y6_tidal_degree200',...
	'y1_loading_degree002','y2_loading_degree002','y3_loading_degree002','y4_loading_degree002','y5_loading_degree002','y6_loading_degree002',...
	'y1_loading_degree020','y2_loading_degree020','y3_loading_degree020','y4_loading_degree020','y5_loading_degree020','y6_loading_degree020',...
	'y1_loading_degree200','y2_loading_degree200','y3_loading_degree200','y4_loading_degree200','y5_loading_degree200','y6_loading_degree200',...
	}; 
field_tolerances={...
    3e-7, 3e-7, 3e-7, 1e-7, 6e-8, 9e-7,...
    2e-7, 7e-8, 3e-7, 9e-8, 9e-10, 8e-10,...
    2e-8, 4e-8, 4e-7, 3e-8, 2e-10, 1e-10,...
    4e-6, 1e-6, 4e-6, 3e-6, 8e-7, 2e-6,...
    2e-6, 1e-7, 5e-6, 3e-7, 2e-7, 2e-7,...
    2e-6, 9e-10, 5e-5, 3e-8, 5e-7, 2e-9...
	}; 
field_values={...
	y1_tidal_degree002,y2_tidal_degree002,y3_tidal_degree002,y4_tidal_degree002,y5_tidal_degree002,y6_tidal_degree002,...
	y1_tidal_degree020,y2_tidal_degree020,y3_tidal_degree020,y4_tidal_degree020,y5_tidal_degree020,y6_tidal_degree020,...
	y1_tidal_degree200,y2_tidal_degree200,y3_tidal_degree200,y4_tidal_degree200,y5_tidal_degree200,y6_tidal_degree200,...
	y1_loading_degree002,y2_loading_degree002,y3_loading_degree002,y4_loading_degree002,y5_loading_degree002,y6_loading_degree002,...
	y1_loading_degree020,y2_loading_degree020,y3_loading_degree020,y4_loading_degree020,y5_loading_degree020,y6_loading_degree020,...
	y1_loading_degree200,y2_loading_degree200,y3_loading_degree200,y4_loading_degree200,y5_loading_degree200,y6_loading_degree200,...
	};

