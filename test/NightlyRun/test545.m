%Test Name: PigTranStochasticforcingCovariance
md=triangle(model(),'../Exp/Pig.exp',10000.);
md=setmask(md,'../Exp/PigShelves.exp','../Exp/PigIslands.exp');
md=parameterize(md,'../Par/Pig.par');
md=setflowequation(md,'SSA','all');
md.timestepping.start_time = 0;
md.timestepping.time_step  = 1;
md.timestepping.final_time = 10;

%Basin separation TF
idb_tf  = zeros(md.mesh.numberofelements,1);
iid1    = find(md.mesh.x<=-1.6e6);
for ii=1:md.mesh.numberofelements
    for vertex=1:3
        if any(iid1==md.mesh.elements(ii,vertex)) %one vertex in basin 1
            idb_tf(ii) = 1;
        end
    end
    if idb_tf(ii)==0 %no vertex was found in basin 1
        idb_tf(ii) = 2;
    end
end
% Basin separation default
idb_df = zeros(md.mesh.numberofelements,1);
iid1   = find(md.mesh.x<=-1.62e6);
for ii=1:md.mesh.numberofelements
    for vertex=1:3
        if any(iid1==md.mesh.elements(ii,vertex)) %one vertex in basin 1
            idb_df(ii) = 1;
        end
    end
    if idb_df(ii)==0 %no vertex was found in basin 1
        idb_df(ii) = 2;
    end
end
% Dimensionalities
nb_tf    = 2;
nb_clv   = 2;
nb_flmlt = 2;

%Calving parameters
md.mask.ice_levelset = 1e4*(md.mask.ice_levelset + 0.5);
md.calving.calvingrate = 0.3*ones(md.mesh.numberofvertices,1);
md.levelset.spclevelset = NaN(md.mesh.numberofvertices,1);
md.levelset.migration_max = 10.0; %avoid fast advance/retreat of the front
%%% Frontal forcing parameters %%%
md.frontalforcings=frontalforcingsrignotarma();
md.frontalforcings.num_basins              = nb_tf;
md.frontalforcings.basin_id                = idb_tf;
% Polynomial params %
md.frontalforcings.num_params        = 1; %only a constant term
md.frontalforcings.num_breaks        = 0; %no breakpoint
constval                             = [2.5;0.5];
md.frontalforcings.polynomialparams  = constval;
% No monthly effects: do nothing %
% ARMA model parameters %
md.frontalforcings.ar_order        = 3;
md.frontalforcings.ma_order        = 2;
md.frontalforcings.arma_timestep   = 2; %timestep of the ARMA model [yr]
md.frontalforcings.arlag_coefs     = [[0.1,-0.1,0.01];[0.2,-0.2,0.1]]; %autoregressive parameters
md.frontalforcings.malag_coefs     = [[0.1,0.0];[0.0,0.1]]; %moving-average parameters
% No ARMA model of subglacial discharge: simply specify values at vertices %
md.frontalforcings.subglacial_discharge = 10*ones(md.mesh.numberofvertices,1);

% Floating Ice Melt parameters
md.basalforcings.floatingice_melting_rate = 0.1*ones(md.mesh.numberofvertices,1);


% Covariance matrix
covtf       = 1e-4*eye(nb_tf);
covclv      = 1e-1*eye(nb_clv);
covclv(1,1) = 1/10*covclv(1,1);
covflmlt    = 0.05*eye(nb_flmlt);
covglob     = blkdiag(covtf,covclv,covflmlt);

%Stochastic forcing
md.stochasticforcing.isstochasticforcing = 1;
md.stochasticforcing.fields              = [{'FrontalForcingsRignotarma'},{'DefaultCalving'},{'FloatingMeltRate'}];
md.stochasticforcing.defaultdimension    = 2;
md.stochasticforcing.default_id          = idb_df;
md.stochasticforcing.covariance          = covglob; %global covariance among- and between-fields
md.stochasticforcing.randomflag          = 0; %determines true/false randomness

md.transient.ismovingfront   = 1;
md.transient.isgroundingline = 1;

md.transient.requested_outputs = {'default','CalvingCalvingrate','CalvingMeltingrate','BasalforcingsFloatingiceMeltingRate'};
md.cluster=generic('name',oshostname(),'np',2);
md=solve(md,'Transient');

%Fields and tolerances to track changes
field_names ={...
   'Vx1' ,'Vy1' ,'Vel1' ,'Thickness1' , 'MaskIceLevelset1', 'CalvingCalvingrate1', 'CalvingMeltingrate1', 'BasalforcingsFloatingiceMeltingRate1',...
   'Vx5' ,'Vy5' ,'Vel5' ,'Thickness5' , 'MaskIceLevelset5', 'CalvingCalvingrate5', 'CalvingMeltingrate5', 'BasalforcingsFloatingiceMeltingRate5',...
   'Vx10' ,'Vy10' ,'Vel10' ,'Thickness10' , 'MaskIceLevelset10', 'CalvingCalvingrate10', 'CalvingMeltingrate10', 'BasalforcingsFloatingiceMeltingRate10',...
   };
field_tolerances={...
   1e-11,2e-11,2e-11,1e-11,1e-9,1e-10,1e-10,1e-10,...
   2e-11,1e-11,1e-11,9e-11,2e-9,1e-10,1e-10,1e-10,...
   2e-6,1e-6,1e-6,1e-6,5e-6,1e-6,1e-6,1e-6,...
   };
field_values={...
   (md.results.TransientSolution(1).Vx),...
   (md.results.TransientSolution(1).Vy),...
   (md.results.TransientSolution(1).Vel),...
   (md.results.TransientSolution(1).Thickness),...
   (md.results.TransientSolution(1).MaskIceLevelset),...
   (md.results.TransientSolution(1).CalvingCalvingrate),...
   (md.results.TransientSolution(1).CalvingMeltingrate),...
   (md.results.TransientSolution(1).BasalforcingsFloatingiceMeltingRate),...
   (md.results.TransientSolution(5).Vx),...
   (md.results.TransientSolution(5).Vy),...
   (md.results.TransientSolution(5).Vel),...
   (md.results.TransientSolution(5).Thickness),...
   (md.results.TransientSolution(5).MaskIceLevelset),...
   (md.results.TransientSolution(5).CalvingCalvingrate),...
   (md.results.TransientSolution(5).CalvingMeltingrate),...
   (md.results.TransientSolution(5).BasalforcingsFloatingiceMeltingRate),...
	(md.results.TransientSolution(10).Vx),...
	(md.results.TransientSolution(10).Vy),...
	(md.results.TransientSolution(10).Vel),...
	(md.results.TransientSolution(10).Thickness),...
	(md.results.TransientSolution(10).MaskIceLevelset),...
	(md.results.TransientSolution(10).CalvingCalvingrate),...
	(md.results.TransientSolution(10).CalvingMeltingrate),...
	(md.results.TransientSolution(10).BasalforcingsFloatingiceMeltingRate),...
	};
