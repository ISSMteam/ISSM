%Test Name: PigTranARMAandStochasticforcings 
md=triangle(model(),'../Exp/Pig.exp',8000.);
md=setmask(md,'../Exp/PigShelves.exp','../Exp/PigIslands.exp');
md=parameterize(md,'../Par/Pig.par');
md=setflowequation(md,'SSA','all');
md.timestepping.start_time = 0;
md.timestepping.time_step  = 1;
md.timestepping.final_time = 10;

%Basin separation
idb     = zeros(md.mesh.numberofelements,1);
iid1    = find(md.mesh.x>=-1.6e6);
for ii=1:md.mesh.numberofelements
    for vertex=1:3
        if any(iid1==md.mesh.elements(ii,vertex)) %one vertex in basin 1
            idb(ii) = 1;
        end
    end
    if idb(ii)==0 %no vertex was found in basin 1
        idb(ii) = 2;
    end
end
nb_bas = 2;

%SMB
numparams               = 1;
numbreaks               = 0;
intercept               = [0.5;0.01];
polynomialparams        = intercept;
datebreaks              = NaN;
md.smb                  = SMBarma();
md.smb.num_basins       = nb_bas; %number of basins
md.smb.basin_id         = idb; %prescribe basin ID number to elements
md.smb.num_params       = numparams; %number of parameters in the polynomial
md.smb.num_breaks       = numbreaks; %number of breakpoints
md.smb.polynomialparams = polynomialparams;
md.smb.datebreaks       = datebreaks;
md.smb.ar_order         = 4;
md.smb.ma_order         = 4;
md.smb.arma_timestep    = 2.0; %timestep of the ARMA model [yr]
md.smb.arlag_coefs      = [[0.2,0.1,0.05,0.01];[0.4,0.2,-0.2,0.1]];
md.smb.malag_coefs      = [[0.1,0.1,0.2,0.3];[0.5,0.8,1.3,2.4]];

%Calving
md.mask.ice_levelset           = 1e4*(md.mask.ice_levelset + 0.5);
md.calving.calvingrate         = 0.1*ones(md.mesh.numberofvertices,1);
md.levelset.spclevelset        = NaN(md.mesh.numberofvertices,1);
md.levelset.migration_max      = 10.0;
md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices,1);

% Basal forcing implementation
numparams                         = 2;
numbreaks                         = 1;
intercept                         = [3.0,4.0;1.0,0.5];
trendlin                          = [0.0,0.1;0,0];
polynomialparams                  = cat(3,intercept,trendlin);
datebreaks                        = [6;7];
md.basalforcings                  = linearbasalforcingsarma();
md.basalforcings.num_basins       = nb_bas; %number of basins
md.basalforcings.basin_id         = idb; %prescribe basin ID number to elements
md.basalforcings.polynomialparams = polynomialparams;
md.basalforcings.datebreaks       = datebreaks;
md.basalforcings.num_params       = numparams; %number of parameters in the polynomial
md.basalforcings.num_breaks       = numbreaks; %number of breakpoints
md.basalforcings.ar_order         = 1;
md.basalforcings.ma_order         = 1;
md.basalforcings.arma_timestep    = 1.0; %timestep of the ARMA model [yr]
md.basalforcings.arlag_coefs      = [0.0;0.1];
md.basalforcings.malag_coefs      = [0.55;0.34];
md.basalforcings.deepwater_elevation       = [-1000,-1520];
md.basalforcings.upperwater_elevation      = [0,-50];
md.basalforcings.upperwater_melting_rate   = [0.0,0.0];
md.basalforcings.groundedice_melting_rate  = zeros(md.mesh.numberofvertices,1);

% Covariance matrix
covsmb      = 3*eye(nb_bas);
covclv      = 1e-1*eye(nb_bas);
covclv(1,1) = 1/10*covclv(1,1);
covdwm      = 400*eye(nb_bas);
covglob     = blkdiag(covsmb,covclv,covdwm);

% Stochastic forcing
md.stochasticforcing.isstochasticforcing = 1;
md.stochasticforcing.fields              = [{'SMBarma'},{'DefaultCalving'},{'BasalforcingsDeepwaterMeltingRatearma'}];
md.stochasticforcing.defaultdimension    = 2;
md.stochasticforcing.default_id          = idb;
md.stochasticforcing.covariance          = covglob; %global covariance among- and between-fields
md.stochasticforcing.randomflag          = 0; %determines true/false randomness

md.transient.ismovingfront     = 1;
md.transient.requested_outputs = {'default','SmbMassBalance','BasalforcingsFloatingiceMeltingRate','BasalforcingsSpatialDeepwaterMeltingRate'};
md.transient.isstressbalance = 1;
md.transient.ismasstransport = 1;
md.transient.issmb           = 1;
md.transient.isthermal       = 0;
md.transient.isgroundingline = 1;

md.cluster=generic('name',oshostname(),'np',2);
md=solve(md,'Transient');

%Fields and tolerances to track changes
field_names ={...
   'Vx1' ,'Vy1' ,'Vel1' ,'Thickness1', 'SmbMassBalance1', 'BasalforcingsFloatingiceMeltingRate1', 'BasalforcingsSpatialDeepwaterMeltingRate1',...
   'Vx2' ,'Vy2' ,'Vel2' ,'Thickness2', 'SmbMassBalance2' ,'BasalforcingsFloatingiceMeltingRate2', 'BasalforcingsSpatialDeepwaterMeltingRate2',...
   'Vx3' ,'Vy3' ,'Vel3' ,'Thickness3', 'SmbMassBalance3' ,'BasalforcingsFloatingiceMeltingRate3', 'BasalforcingsSpatialDeepwaterMeltingRate3',...
   };
field_tolerances={...
   1e-11,1e-11,2e-11,1e-11,1e-10,1e-9,1e-10,...
   1e-11,1e-11,2e-11,9e-11,1e-10,1e-9,1e-10,...
   2e-10,2e-10,2e-10,1e-10,1e-10,1e-9,1e-10,...
   };
field_values={...
   (md.results.TransientSolution(1).Vx),...
   (md.results.TransientSolution(1).Vy),...
   (md.results.TransientSolution(1).Vel),...
   (md.results.TransientSolution(1).Thickness),...
   (md.results.TransientSolution(1).SmbMassBalance),...
   (md.results.TransientSolution(1).BasalforcingsFloatingiceMeltingRate),...
   (md.results.TransientSolution(1).BasalforcingsSpatialDeepwaterMeltingRate),...
   (md.results.TransientSolution(5).Vx),...
   (md.results.TransientSolution(5).Vy),...
   (md.results.TransientSolution(5).Vel),...
   (md.results.TransientSolution(5).Thickness),...
   (md.results.TransientSolution(5).SmbMassBalance),...
   (md.results.TransientSolution(5).BasalforcingsFloatingiceMeltingRate),...
   (md.results.TransientSolution(5).BasalforcingsSpatialDeepwaterMeltingRate),...
	(md.results.TransientSolution(10).Vx),...
	(md.results.TransientSolution(10).Vy),...
	(md.results.TransientSolution(10).Vel),...
	(md.results.TransientSolution(10).Thickness),...
   (md.results.TransientSolution(10).SmbMassBalance),...
	(md.results.TransientSolution(10).BasalforcingsFloatingiceMeltingRate),...
   (md.results.TransientSolution(10).BasalforcingsSpatialDeepwaterMeltingRate),...
	};
