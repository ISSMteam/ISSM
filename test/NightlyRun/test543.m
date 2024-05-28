%Test Name: PigTranRignotarma
md=triangle(model(),'../Exp/Pig.exp',10000.);
md=setmask(md,'../Exp/PigShelves.exp','../Exp/PigIslands.exp');
md=parameterize(md,'../Par/Pig.par');
md=setflowequation(md,'SSA','all');
md.timestepping.start_time = 0;
md.timestepping.time_step  = 0.05;
md.timestepping.final_time = 2;

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
% Dimensionalities
nb_tf    = 2;

%Calving parameters
md.mask.ice_levelset = 1e4*(md.mask.ice_levelset + 0.5);
md.calving.calvingrate = 0*ones(md.mesh.numberofvertices,1);
md.levelset.spclevelset = NaN(md.mesh.numberofvertices,1);
md.levelset.migration_max = 10.0; %avoid fast advance/retreat of the front
%%% Frontal forcing parameters %%%
md.frontalforcings=frontalforcingsrignotarma();
% Polynomial params %
numparams        = 3;
numbreaks        = 2;
intercept        = [2.5,2.0,0.1;0.5,0.5,1.5];
trendlin         = [-0.5,-0.2,0.1;0,0,0];
trendquad        = [0,0.0,0;0.1,0.1,0.1];
datebreaks       = [0.5,1.5;0.5,1.5];
polynomialparams = cat(numparams,intercept,trendlin,trendquad);
% Monthly effects params %
numbreaksM       = 1;
intcpsMp0        = [-0.5,-0.5,0,0,0,0,0.5,0.5,0,0,0,0;
                    -1.0,-1.0,0,0,0,0,1.0,1.0,0,0,0,0];
intcpsMp1        = [-0.25,-0.25,0,0,0,0,0.,0.,0,0,0,0;
                    -0.1,-0.1,0,0,0,0,0.1,0.1,0,0,0,0];
intcpsM          = cat(3,intcpsMp0,intcpsMp1);
trendsMp0        = [0,0,0,0,0,0,0.,0.0,0,0,0,0;
                    0.0,0.0,0,-0.0,0,0,0.0,0.0,0,0,0,0];
trendsMp1        = [0,-0.12,0,0,0,0,0.,0.0,0,0.0,0,0;
                    0.0,-0.1,0,-0.0,0,0,0.0,0.0,0,0,0,0];
trendsM          = cat(3,trendsMp0,trendsMp1);
datebreaksM      = [1;1];
% Subglacial discharge parameters %
isdischargearma            = 1;
sd_ar_order                = 1;
sd_ma_order                = 1;
sd_num_breaks              = 1;
sd_num_params              = 2;
sd_arma_timestep           = 1;
sd_arlag_coefs             = [0.95;0.95];
sd_malag_coefs             = [0;0];
sd_datebreaks              = [1;1];
sd_monthlyfrac             = [0,0,0,0,0,0,0.5,0.5,0,0,0,0;
                              0,0,0,0,0,0,0.5,0.5,0,0,0,0];
sd_const                   = [50000,70000.0;8000,10000.0];
sd_trend                   = [0,10000;0,0];
sd_polyparam               = cat(3,sd_const,sd_trend);

md.frontalforcings.num_basins              = nb_tf;
md.frontalforcings.basin_id                = idb_tf;
md.frontalforcings.num_params              = numparams; %number of parameters in the polynomial
md.frontalforcings.num_breaks              = numbreaks; %number of breakpoints
md.frontalforcings.subglacial_discharge    = 0.01*ones(md.mesh.numberofvertices,1);
md.frontalforcings.polynomialparams        = polynomialparams;
md.frontalforcings.datebreaks              = datebreaks;
md.frontalforcings.ar_order                = 4;
md.frontalforcings.ma_order                = 2;
md.frontalforcings.arma_timestep           = 2; %timestep of the ARMA model [yr]
md.frontalforcings.arlag_coefs             = [[0.1,-0.1,0.01,-0.01];[0.2,-0.2,0.1,0.0]]; %autoregressive parameters
md.frontalforcings.malag_coefs             = [[0.1,0.0];[0.0,0.1]]; %moving-average parameters
md.frontalforcings.monthlyvals_numbreaks   = numbreaksM;
md.frontalforcings.monthlyvals_intercepts  = intcpsM;
md.frontalforcings.monthlyvals_trends      = trendsM;
md.frontalforcings.monthlyvals_datebreaks  = datebreaksM;
md.frontalforcings.isdischargearma         = isdischargearma;
if(isdischargearma==0)
	md.frontalforcings.subglacial_discharge    = 0.01*ones(md.mesh.numberofvertices,1);
else
    md.frontalforcings.sd_num_breaks         = sd_num_breaks;
    md.frontalforcings.sd_num_params         = sd_num_params;
    md.frontalforcings.sd_ar_order           = sd_ar_order;
    md.frontalforcings.sd_ma_order           = sd_ma_order;
    md.frontalforcings.sd_arma_timestep      = sd_arma_timestep;
    md.frontalforcings.sd_arlag_coefs        = sd_arlag_coefs;
    md.frontalforcings.sd_malag_coefs        = sd_malag_coefs;
    md.frontalforcings.sd_datebreaks         = sd_datebreaks;
    md.frontalforcings.sd_monthlyfrac        = sd_monthlyfrac;
    md.frontalforcings.sd_polynomialparams   = sd_polyparam;
end
% Floating Ice Melt parameters
md.basalforcings.floatingice_melting_rate = 0*ones(md.mesh.numberofvertices,1);


% Covariance matrix
covtf       = 1e-4*eye(nb_tf);
covsd       = 1e3*eye(nb_tf);
covglob     = blkdiag(covtf,covsd);

%Stochastic forcing
md.stochasticforcing.isstochasticforcing = 1;
md.stochasticforcing.fields              = [{'FrontalForcingsRignotarma'},{'FrontalForcingsSubglacialDischargearma'}];
md.stochasticforcing.defaultdimension    = 2;
md.stochasticforcing.default_id          = idb_tf;
md.stochasticforcing.covariance          = covglob; %global covariance among- and between-fields
md.stochasticforcing.randomflag          = 0; %determines true/false randomness

md.transient.ismovingfront   = 1;
md.transient.isgroundingline = 1;

md.transient.requested_outputs = {'default','CalvingMeltingrate'};
md.cluster=generic('name',oshostname(),'np',2);
md=solve(md,'Transient');

%Fields and tolerances to track changes
field_names ={...
   'Vx1' ,'Vy1' ,'Vel1' ,'Thickness1' ,'MaskIceLevelset1' ,'CalvingMeltingrate1' ,...
   'Vx2' ,'Vy2' ,'Vel2' ,'Thickness2' ,'MaskIceLevelset2' ,'CalvingMeltingrate2' ,...
   'Vx10','Vy10','Vel10','Thickness10','MaskIceLevelset10','CalvingMeltingrate10',...
   };
field_tolerances={...
   1e-11,2e-11,2e-11,1e-11,1e-9,1e-10,...
   2e-11,1e-11,1e-11,9e-11,2e-9,1e-10,...
   2e-6,1e-6,1e-6,1e-6,5e-6,1e-6,...
   };
field_values={...
   (md.results.TransientSolution(1).Vx),...
   (md.results.TransientSolution(1).Vy),...
   (md.results.TransientSolution(1).Vel),...
   (md.results.TransientSolution(1).Thickness),...
   (md.results.TransientSolution(1).MaskIceLevelset),...
   (md.results.TransientSolution(1).CalvingMeltingrate),...
   (md.results.TransientSolution(20).Vx),...
   (md.results.TransientSolution(20).Vy),...
   (md.results.TransientSolution(20).Vel),...
   (md.results.TransientSolution(20).Thickness),...
   (md.results.TransientSolution(20).MaskIceLevelset),...
   (md.results.TransientSolution(20).CalvingMeltingrate),...
	(md.results.TransientSolution(40).Vx),...
	(md.results.TransientSolution(40).Vy),...
	(md.results.TransientSolution(40).Vel),...
	(md.results.TransientSolution(40).Thickness),...
	(md.results.TransientSolution(40).MaskIceLevelset),...
	(md.results.TransientSolution(40).CalvingMeltingrate),...
	};
