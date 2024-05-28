%Test Name: SHslc 
% spherical-harmonic (SH) approach for solving the sea-level equation. 
% reference: Adhikari et al., 2019, ESSD, https://doi.org/10.5194/essd-11-629-2019 

% parameters {{{ 
	% MEALpix parameters 
	num  = 5;   % choose one of [4-7] :: WARNING :: num=7 takes a while to solve!
	lMax = 60;  % SH of maximum degree
	nSide = 2^num;
	nPix  = 12*nSide^2; % number of pixels 
	num_lm = (1+lMax)^2; % number of SH coefficients 

	% relative tolerance 
	para.rel_tol = 1.0e-5; 

	% love numbers 
	love=lovenumbers('maxdeg',10000,'referenceframe','CM');
	para.loveK = love.k; 
	para.loveH = love.h; 
	para.k2 = love.tk(3);		% tidal Love number 
	para.h2 = love.th(3);		% tidal Love number 
	para.ks = love.tk2secular; 

	md=model();
	para.A = md.solidearth.rotational.equatorialmoi;	% mean equatorial moment of inertia, kg m^2 
	para.C = md.solidearth.rotational.polarmoi;			% polar moment of inertia, kg m^2 
	para.Omega =  md.solidearth.rotational.angularvelocity; % mean rotational velocity of Earth, 1/s 
	para.earth_radius = md.solidearth.planetradius;	% mean radius of the Earth, m 
	para.rho_ocean = md.materials.rho_water;  % ocean water density, kg/m^3
	para.rho_earth = md.materials.earth_density;  % average density of the Earth, kg/m^3
	para.g = md.love.g0; % mean surface gravitational acceleration, m/s^2 
	clearvars md; 

% }}} 

% Compute real spherical harmonics for chosen mealpix resolution {{{
   disp([' === Computing "real" spherical harmonics for chosen resolution ==== ']); 
   clearvars sh  % if there is any! 
   q = 0; 
	pix_center_gmtdir=pix2ang(nSide);
	pix_center=cell2mat(pix_center_gmtdir);
	lat_p=pix_center(1,:)';
	lon_p=pix_center(2,:)';
	for l=0:lMax
		plm = legendre(l,cos(lat_p),'norm'); %nromalized Plm (see Matlab documents) 
		for m=-l:l 
			sh(:,1+q)=plm(abs(m)+1,:)'.*(cos(abs(m).*lon_p)*(m>=0)+sin(abs(m).*lon_p)*(m<0)).*sqrt((2-(m==0))*2); %4-pi
			q=q+1; 
		end 
      disp(['     SH computations : degree ', num2str(l), ' (of ', num2str(lMax), ') completed.' ]);  
	end 
   disp(['     completed! ']); 
	
	% lat,lon for plotting 
	lat_deg=lat_p'*180/pi; 
	lon_deg=lon_p'*180/pi; 
	lat_deg = 90-lat_deg;
	lon_deg(lon_deg>180) = lon_deg(lon_deg>180)-360; 
%}}}

% Compute SH coefficients for ocean function {{{ 
	disp([' === Computing SH coefficients for ocean function ================== ']); 
	ocean = gmtmask(lat_deg',lon_deg')'; 
	land = 1-ocean;	% continental function 
%	oce_lm = funlm(lMax,nPix,ocean,sh);  % compute SH coefficients for ocean function  
	disp(['     completed! ']); 
%}}}

% forcing (weh) field {{{ 
	load('../Data/GRACE_JPL_April2002_WEH.mat');
	lat = repmat(lat',720,1); 
	lon = repmat(lon,1,360); 

	F = scatteredInterpolant(lat(:),lon(:),weh(:)); 
	force = F(lat_deg,lon_deg);
	force(isnan(force))=0; 

	force = force.*land; % only keep loads on land.  

% }}} 

% solve. 
disp([' === Solving WOUT rotational feedbacks =========']); 
rotational_feedback=0; % without rotational feedbacks. 
[bary_0,rsl_0,geoid_0,bed_0] = SHslr(sh,para,force,ocean,rotational_feedback); 

disp([' === Solving WITH rotational feedbacks =========']); 
rotational_feedback=1; % with rotational feedbacks. 
[bary_1,rsl_1,geoid_1,bed_1] = SHslr(sh,para,force,ocean,rotational_feedback); 

%Fields and tolerances to track changes
field_names={'Barystatic_no_rotation','rsl_no_rotation','geoid_no_rotation','bed_no_rotation'...
	'Barystatic_with_rotation','rsl_with_rotation','geoid_with_rotation','bed_with_rotation'};
field_tolerances={3e-7,1e-13,1e-13,1e-13,...
	2e-7,1e-13,1e-13,1e-13};
field_values={bary_0,rsl_0,geoid_0,bed_0,...
	bary_1,rsl_1,geoid_1,bed_1}; 

return; 

% plot. 
% plot {{{ 
	field = (rsl_1-rsl_0)'; 

	res = 1.0; 
	[lat_grid, lon_grid] = meshgrid(linspace(-90,90,180/res), linspace(-180,180,360/res));
	sol_grid = zeros(size(lat_grid)); 

	F = scatteredInterpolant(lat_deg',lon_deg',field); 
	F.Method = 'linear';
	F.ExtrapolationMethod = 'linear'; 

	sol_grid = F(lat_grid, lon_grid);
	sol_grid(isnan(sol_grid))=0; 

	set(0,'DefaultAxesFontSize',18,'DefaultAxesLineWidth',1,'DefaultTextFontSize',18,'DefaultLineMarkerSize',8)
	figure1=figure('Position', [100, 100, 1000, 500]); 
	gcf; load coast; cla; 
	pcolor(lon_grid,lat_grid,sol_grid); shading flat; hold on;
	plot(long,lat,'k'); hold off; 
	c1=colorbar;
	colormap('haxby');
% }}} 


