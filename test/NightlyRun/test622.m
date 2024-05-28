%Test Name:79NorthHydrologyArmapw 
md=triangle(model(),'../Exp/79North.exp',6000.);
md=setmask(md,'../Exp/79NorthShelf.exp','');
md=parameterize(md,'../Par/79North.par');
md=setflowequation(md,'SSA','all');

%Default friction
md.friction             = friction();
md.friction.coefficient = 30*ones(md.mesh.numberofvertices,1);
md.friction.p           = 1*ones(md.mesh.numberofelements,1);
md.friction.q           = 1*ones(md.mesh.numberofelements,1);

% Basin separation default
idb_df = zeros(md.mesh.numberofelements,1);
iid1   = find(md.mesh.y<=-1.08e6);
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
%%% Covariance matrix %%%
covPw  = [1e9,0;0,1e9];
covSMB = [0.1,0;0,0.1];
covGlob = blkdiag(covPw,covSMB);

%%% Hydrology scheme %%%
md.hydrology             = hydrologyarmapw();
md.hydrology.num_basins            = 2;
md.hydrology.basin_id              = idb_df;
md.hydrology.monthlyfactors        = 1*ones(md.hydrology.num_basins,12);
md.hydrology.monthlyfactors(:,1:3) = 0;
md.hydrology.num_params            = 2; %number of parameters in the polynomial
md.hydrology.num_breaks            = 2; %number of breakpoints
termconst                          = [0.5*1e6,0.1*1e6,0.5e6;
                                      0.5*1e6,0.1*1e6,0.5e6];
termtrend                          = [1*1e5,0,0.;
                                      0,1*1e5,0];
md.hydrology.polynomialparams      = cat(3,termconst,termtrend);
md.hydrology.datebreaks            = [20,40;20,40];
md.hydrology.arma_timestep         = 1;
md.hydrology.ar_order              = 1;
md.hydrology.ma_order              = 1;
md.hydrology.arlag_coefs           = [0.98;0.98];
md.hydrology.malag_coefs           = [0;0];

% SMB
md.smb = SMBarma();
md.smb.num_basins            = 2;
md.smb.basin_id              = idb_df;
md.smb.num_breaks            = 0;
md.smb.num_params            = 1;
md.smb.polynomialparams      = 0*[0.5;0.2];
md.smb.ar_order              = 1;
md.smb.ma_order              = 1;
md.smb.arlag_coefs           = [0;0];
md.smb.malag_coefs           = [0;0];
md.smb.arma_timestep         = 1.0;

%%% Stochastic forcing %%%
md.stochasticforcing.isstochasticforcing = 1;
md.stochasticforcing.fields              = [{'FrictionWaterPressure'},{'SMBarma'}];
md.stochasticforcing.defaultdimension    = 2;
md.stochasticforcing.default_id          = idb_df;
md.stochasticforcing.covariance          = covGlob; %global covariance
md.stochasticforcing.stochastictimestep  = 1; %time step of stochastic forcing
md.stochasticforcing.randomflag          = 0; %determines true/false randomness

md.transient.issmb              = 1;
md.transient.ismasstransport    = 1;
md.transient.isstressbalance    = 1;
md.transient.isthermal          = 0;
md.transient.isgroundingline    = 0;
md.transient.ishydrology        = 1;

md.transient.requested_outputs = {'default','SmbMassBalance','FrictionWaterPressure'};
md.timestepping.start_time = 0;
md.timestepping.time_step  = 1.0/12;
md.timestepping.final_time = 2;
md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Transient');

%Fields and tolerances to track changes
field_names      = {'Vel1','Thickness1','SmbMassBalance1','FrictionWaterPressure1',...
                    'Vel12','Thickness12','SmbMassBalance12','FrictionWaterPressure12',...
                    'Vel24','Thickness24','SmbMassBalance24','FrictionWaterPressure24'};
field_tolerances={2e-10,2e-10,2e-10,2e-10,...
                  4e-10,4e-10,4e-10,4e-10,...
                  8e-10,8e-10,8e-10,8e-10};
              
field_values={...
    (md.results.TransientSolution(1).Vel),...
    (md.results.TransientSolution(1).Thickness),...
    (md.results.TransientSolution(1).SmbMassBalance),...
    (md.results.TransientSolution(1).FrictionWaterPressure),...
    (md.results.TransientSolution(12).Vel),...
    (md.results.TransientSolution(12).Thickness),...
    (md.results.TransientSolution(12).SmbMassBalance),...
    (md.results.TransientSolution(12).FrictionWaterPressure),...
    (md.results.TransientSolution(24).Vel),...
    (md.results.TransientSolution(24).Thickness),...
    (md.results.TransientSolution(24).SmbMassBalance),...
    (md.results.TransientSolution(24).FrictionWaterPressure),...
    };



