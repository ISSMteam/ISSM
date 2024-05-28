%Test Name: 79NorthStochFrictionWaterPressure
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
%Covariance matrix
covPw      = 0.5e10*eye(2);
covPw(1,1) = 1.5*covPw(1,1);

%Stochastic forcing
md.stochasticforcing.isstochasticforcing = 1;
md.stochasticforcing.fields              = [{'FrictionWaterPressure'}];
md.stochasticforcing.defaultdimension    = 2;
md.stochasticforcing.default_id          = idb_df;
md.stochasticforcing.covariance          = covPw; %global covariance
md.stochasticforcing.randomflag          = 0; %determines true/false randomness

md.transient.issmb              = 0;
md.transient.ismasstransport    = 1;
md.transient.isstressbalance    = 1;
md.transient.isthermal          = 0;
md.transient.isgroundingline    = 0;

md.transient.requested_outputs = {'default','FrictionWaterPressure'};
md.timestepping.start_time = 0;
md.timestepping.time_step  = 1;
md.timestepping.final_time = 5;
md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Transient');


%Fields and tolerances to track changes
field_names      = {'Vx1','Vy1','Vel1','Thickness1','FrictionWaterPressure1',...
                    'Vx2','Vy2','Vel2','Thickness2','FrictionWaterPressure2',...
                    'Vx5','Vy5','Vel5','Thickness5','FrictionWaterPressure5'};
field_tolerances={2e-10,2e-10,2e-10,2e-10,2e-10,...
                  4e-10,4e-10,4e-10,4e-10,4e-10,...
                  8e-10,8e-10,8e-10,8e-10,8e-10};
              
field_values={...
	(md.results.TransientSolution(1).Vx),...
    (md.results.TransientSolution(1).Vy),...
    (md.results.TransientSolution(1).Vel),...
    (md.results.TransientSolution(1).Thickness),...
    (md.results.TransientSolution(1).FrictionWaterPressure),...
    (md.results.TransientSolution(2).Vx),...
    (md.results.TransientSolution(2).Vy),...
    (md.results.TransientSolution(2).Vel),...
    (md.results.TransientSolution(2).Thickness),...
    (md.results.TransientSolution(2).FrictionWaterPressure),...
    (md.results.TransientSolution(5).Vx),...
    (md.results.TransientSolution(5).Vy),...
    (md.results.TransientSolution(5).Vel),...
    (md.results.TransientSolution(5).Thickness),...
    (md.results.TransientSolution(5).FrictionWaterPressure),...
    };

