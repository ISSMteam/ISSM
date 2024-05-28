
%%%
% Tutorial for StISSM
%%%

steps = [1];

if any(steps==1)
   % Mesh parameters 
   domain    = ['./DomainOutline.exp'];
   hinit     = 5000; % element size for the initial mesh
   hmax      = 40000; % maximum element size of the final mesh
   hmin      = 4000; % minimum element size of the final mesh
   gradation = 1.7; % maximum size ratio between two neighboring elements
   err       = 8; % maximum error between interpolated and control field

   % Generate an initial uniform mesh (resolution = hinit m)
   md = bamg(model,'domain',domain,'hmax',hinit);

   % Get necessary data to build up the velocity grid
   nsidc_vel  = '../Data/Antarctica_ice_velocity.nc';
   xmin       = strsplit(ncreadatt(nsidc_vel,'/','xmin')); 
   xmin       = str2num(xmin{2});
   ymax       = strsplit(ncreadatt(nsidc_vel,'/','ymax'));
   ymax       = str2num(ymax{2});
   spacing    = strsplit(ncreadatt(nsidc_vel,'/','spacing')); 
   spacing    = str2num(spacing{2});
   nx         = double(ncreadatt(nsidc_vel,'/','nx'));
   ny         = double(ncreadatt(nsidc_vel,'/','ny'));
   vx         = double(ncread(nsidc_vel,'vx'));
   vy         = double(ncread(nsidc_vel,'vy'));

   % Build the coordinates
   x = xmin+(0:1:nx)'*spacing;
   y = (ymax-ny*spacing)+(0:1:ny)'*spacing;

   % Interpolate velocities onto coarse mesh
   vx_obs   = InterpFromGridToMesh(x,y,flipud(vx'),md.mesh.x,md.mesh.y,0);
   vy_obs   = InterpFromGridToMesh(x,y,flipud(vy'),md.mesh.x,md.mesh.y,0);
   vel_obs  = sqrt(vx_obs.^2+vy_obs.^2);
   clear vx vy x y;

   % Adapt the mesh to minimize error in velocity interpolation
   md = bamg(md,'hmax',hmax,'hmin',hmin,'gradation',gradation,'field',vel_obs,'err',err);

   % Plot and save model
   plotmodel(md,'data','mesh')
   save ./Models/PIG_StISSM_Mesh_generation md;
   
end

if any(steps==2) %Masks #2 
   md = loadmodel('./Models/PIG_StISSM_Mesh_generation');

   % Load SeaRISe dataset for Antarctica  http://websrv.cs.umt.edu/isis/index.php/Present_Day_Antarctica
   searise = '../Data/Antarctica_5km_withshelves_v0.75.nc';

   % Read thickness mask from SeaRISE
   x1      = double(ncread(searise,'x1'));
   y1      = double(ncread(searise,'y1'));
   thkmask = double(ncread(searise,'thkmask'));

   % Interpolate onto our mesh vertices
   groundedice                 = double(InterpFromGridToMesh(x1,y1,thkmask',md.mesh.x,md.mesh.y,0));
   groundedice(groundedice<=0) = -1;
   clear thkmask;

   % Fill in the md.mask structure
   md.mask.ocean_levelset = groundedice; %ice is grounded for mask equal one
   md.mask.ice_levelset   = -1*ones(md.mesh.numberofvertices,1); %ice is present when negatvie

   plotmodel(md,'data',md.mask.ocean_levelset,'title','grounded/floating','data',md.mask.ice_levelset,'title','ice/no-ice')

   save ./Models/PIG_StISSM_SetMask md;
end

if any(steps==3) %Parameterization #3 
   
   md = loadmodel('./Models/PIG_StISSM_SetMask');
   md = setflowequation(md,'SSA','all');
   md = parameterize(md,'./PigStISSM.par');

   save ./Models/PIG_StISSM_Parameterization md;
end

if any(steps==4) %Stochastic SMB
    md = loadmodel('./Models/PIG_StISSM_Parameterization');

    % Create the different subdomains for SMB %
    ymax  = max(md.mesh.y);
    ymin  = min(md.mesh.y);
    xmax  = max(md.mesh.x);
    xmin  = min(md.mesh.x);
    idsmb = zeros(md.mesh.numberofelements,1); %subdomain ID is defined for each element
    iid1  = find(md.mesh.x>=(xmax-2/3*(xmax-xmin))); %vertices in subdomain 1
    iid2  = find(md.mesh.x<(xmax-2/3*(xmax-xmin))); %vertices in subdomain 2
    for ii=1:md.mesh.numberofelements
        for vertex=1:3
            if any(iid1==md.mesh.elements(ii,vertex)) %one vertex in subdomain 1
                idsmb(ii) = 1;
            end
        end
        if idsmb(ii)==0 %no vertex was found in subdomain 1
            idsmb(ii) = 2;
        end
    end
    
    % SMBarma implementation 
    md.smb = SMBarma();
    md.smb.num_basins = 2; %we use two different subdomains
    md.smb.basin_id   = idsmb; %element subdomain IDs
    md.smb.num_breaks = 1; %1 breakpoint in the piecewise polynomial
    md.smb.datebreaks = [5;5]; %breakpoint occurs at year 5 in both subdomains
    md.smb.num_params = 2; %use a constant and a linear trend for the piecwise polynomial    
    constsmb    = [0.5,0.2;0.3,0.5]; %constant SMB term for subdomains pre- and post-breakpoint [m yr^-1] 
    trendsmb    = [0,0;0.01,0.001]; %trend is SMB for subdomains pre- and post-breakpoint [m yr^-2] 
    md.smb.polynomialparams = cat(3,constsmb,trendsmb); %concatenating const and trend along a 3rd dimension
    md.smb.ar_order         = 1; %first-order AR
    md.smb.ma_order         = 1; %first-order MA
    md.smb.arlag_coefs      = [0.3;0]; %AR coefficients in each subdomain
    md.smb.malag_coefs      = zeros(md.smb.num_basins,md.smb.ma_order); %all zeros is equivalent to MA order 0
    md.smb.arma_timestep    = 1.0; %yearly ARMA  
    md.smb.elevationbins    = [300,1000;300,1000]; %elevations separating different lapse rate values [m]
    md.smb.lapserates       = 1e-2*[0.03,0.01,-1e-4;0.02,0.02,-1e-5]; %lapse rate values [m ice eq. m^-1 yr^-1]
    md.smb.refelevation     = [500,500]; %elevations at which the SMBarma calculated values apply (i.e., before using lapse rates)
    
    % Set-up the covariance matrix
    sdevSMB1         = 0.01; %low standard deviation in subdomain 1 [m ice eq. yr^-1]
    sdevSMB2         = 0.2; %higher variability in subdomain 2 [m ice eq. yr^-1]
    correlationSMB   = [1.0,0.5;0.5,1.0]; %moderate correlation between the subdomains
    covarianceSMB    = diag([sdevSMB1,sdevSMB2])*correlationSMB*diag([sdevSMB1,sdevSMB2]); %covariance matrix [(m ice eq. yr^-1)^2]
    
    % Stochasticforcing implementation
    md.stochasticforcing.isstochasticforcing = 1; %activate stochasticity
    md.stochasticforcing.fields              = [{'SMBarma'}];
    md.stochasticforcing.covariance          = covarianceSMB; %prescribe the SMB covariance
    md.stochasticforcing.stochastictimestep  = 1.0; %yearly stochasticity
    
    save ./Models/PIG_StISSM_StochSMB md;
    
end

if any(steps==5) % Transient Run #1 

   md = loadmodel('./Models/PIG_StISSM_StochSMB');

   md.inversion.iscontrol         = 0;
   md.transient.ismasstransport   = 1;
   md.transient.isstressbalance   = 1;
   md.transient.isgroundingline   = 1;
   md.transient.ismovingfront     = 0;
   md.transient.isthermal         = 0;
   md.verbose.solution            = 1;
   md.timestepping.time_step      = 0.1;
   md.timestepping.final_time     = 10;
   md.transient.requested_outputs = {'default','SmbMassBalance','MaskIceLevelset'};

   md        = solve(md,'Transient');
   timing    = cell2mat({md.results.TransientSolution(:).time});
   fullsmb   = cell2mat({md.results.TransientSolution(:).SmbMassBalance});
   timemeansmb = mean(fullsmb,2);
   timesdevsmb = std(fullsmb,0,2);
   maxsmb      = max(fullsmb,[],'all');
   minsmb      = min(fullsmb,[],'all');
   
   
   plotmodel(md,'figure',1,...
       'data',md.results.TransientSolution(1).SmbMassBalance,'title','SMB at timestep 1','caxis#1',[minsmb,maxsmb],...
       'data',md.results.TransientSolution(end).SmbMassBalance,'title','SMB at last timestep','caxis#2',[minsmb,maxsmb])
   plotmodel(md,'figure',2,...
       'data',timemeansmb,'title','SMB mean over simulation',...
       'data',timesdevsmb,'title','SMB standard deviation over simulation')
       
   seriesSMB1  = fullsmb(1,:);
   seriesSMB71 = fullsmb(71,:);
   figure(3);
   plot(timing,seriesSMB1,'b');
   hold on
   plot(timing,seriesSMB71,'r');
   title('SMB at two different vertices');
   hold off
   
   save ./Models/PIG_StISSM_Transient1 md;
   
end

if any(steps==6) % Set up stochastic calving also
    
    md = loadmodel('./Models/PIG_StISSM_Transient1');

    md.calving.calvingrate         = 20*ones(md.mesh.numberofvertices,1);
    md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices,1);
    md.levelset.spclevelset        = NaN(md.mesh.numberofvertices,1);
    md.levelset.migration_max      = 100.0; %avoid too fast advance/retreat of the front
    md.transient.ismovingfront     = 1;
    
    % Parameterize stochastic calving
    idbasinsCalving = md.smb.basin_id; %same subdomains as for SMB
    sdevClv1         = 0.5; %low standard deviation in subdomain 1 [m yr^-1]
    sdevClv2         = 5; %higher variability in subdomain 2 [m yr^-1]
    correlationClv   = [1.0,0.0;0.0,1.0]; %no correlation between the subdomains
    covarianceClv    = diag([sdevClv1,sdevClv2])*correlationClv*diag([sdevClv1,sdevClv2]); %covariance matrix [(m yr^-1)^2]
    
    % Adjust stochastic forcing class
    md.stochasticforcing.fields     = [{'SMBarma'},{'DefaultCalving'}];
    oldcovarianceSMB                = md.stochasticforcing.covariance;
    covarianceGlobal                = blkdiag(oldcovarianceSMB,covarianceClv); %independence between SMB and calving
    md.stochasticforcing.covariance = covarianceGlobal;
    md.stochasticforcing.defaultdimension = md.smb.num_basins; %2 basins for stochastic default (used for calving, same as SMB)
    md.stochasticforcing.default_id       = idbasinsCalving;
    
    save ./Models/PIG_StISSM_StochCalving md;  
       
end

if any(steps==7) % Transient Run #2

   md = loadmodel('./Models/PIG_StISSM_StochCalving');

   md.verbose.solution            = 1;
   md.timestepping.start_time     = md.timestepping.final_time;
   md.timestepping.time_step      = 0.1;
   md.timestepping.final_time     = md.timestepping.start_time+5;
   md.transient.requested_outputs = {'default','SmbMassBalance','CalvingCalvingrate'};

   % Set up initial conditions from previous transient results
   md.geometry.thickness        = md.results.TransientSolution(end).Thickness;
   md.initialization.vx         = md.results.TransientSolution(end).Vx;
   md.initialization.vy         = md.results.TransientSolution(end).Vy;
   md.initialization.vel        = md.results.TransientSolution(end).Vel;
   md.mask.ocean_levelset       = md.results.TransientSolution(end).MaskOceanLevelset;
   md.initialization.pressure   = md.results.TransientSolution(end).Pressure;
   md.geometry.base             = md.results.TransientSolution(end).Base;
   md.geometry.surface          = md.results.TransientSolution(end).Surface;
   md.mask.ice_levelset         = md.results.TransientSolution(end).MaskIceLevelset;   
   
   md        = solve(md,'Transient');
   md        = loadresultsfromcluster(md);
   timing    = cell2mat({md.results.TransientSolution(:).time});
   fullsmb   = cell2mat({md.results.TransientSolution(:).SmbMassBalance});
   fullclv   = cell2mat({md.results.TransientSolution(:).CalvingCalvingrate});

   timemeansmb = mean(fullsmb,2);
   timesdevsmb = std(fullsmb,0,2);
   maxsmb      = max(fullsmb,[],'all');
   minsmb      = min(fullsmb,[],'all');
   timemeanclv = mean(fullclv,2);
   timesdevclv = std(fullclv,0,2);
   maxclv      = max(fullclv,[],'all');
   minclv      = min(fullclv,[],'all');
   
   plotmodel(md,'figure',1,...
       'data',md.results.TransientSolution(1).SmbMassBalance,'title','SMB at timestep 1','caxis#1',[minsmb,maxsmb],...
       'data',md.results.TransientSolution(end).SmbMassBalance,'title','SMB at last timestep','caxis#2',[minsmb,maxsmb])
   plotmodel(md,'figure',2,...
       'data',timemeansmb,'title','SMB mean over simulation',...
       'data',timesdevsmb,'title','SMB standard deviation over simulation')
   plotmodel(md,'figure',3,...
       'data',timemeanclv,'title','Calving mean over simulation',...
       'data',timesdevclv,'title','Calving standard deviation over simulation')  
   
   seriesSMB1  = fullsmb(1,:);
   seriesSMB71 = fullsmb(71,:);
   figure(4);
   plot(timing,seriesSMB1,'b');
   hold on
   plot(timing,seriesSMB71,'r');
   title('SMB at two different vertices');
   hold off 
   
   seriesclv1  = fullclv(1,:);
   seriesclv71 = fullclv(71,:);
   figure(5);
   plot(timing,seriesclv1,'b');
   hold on
   plot(timing,seriesclv71,'r');
   title('Calving at two different vertices');
   hold off 
   
end








