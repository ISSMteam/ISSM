%Test Name: SquareShelfSMBarma
md=triangle(model(),'../Exp/Square.exp',80000.);
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareShelfConstrained.par');
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',3);
md.transient.requested_outputs={'default','IceVolume','SmbMassBalance'};

ymax = max(md.mesh.y);
xmax = max(md.mesh.x);
%Generate basin IDs for 3 basins
idbasin      = zeros(md.mesh.numberofelements,1);
iid1         = find(md.mesh.y>=2/3*ymax);
iid2         = intersect(find(md.mesh.y<2/3*ymax),find(md.mesh.x>=1/3*xmax));
iid3         = intersect(find(md.mesh.y<2/3*ymax),find(md.mesh.x<1/3*xmax));
for ii=1:md.mesh.numberofelements
    for vertex=1:3
        if any(iid1==md.mesh.elements(ii,vertex)) %one vertex in basin 1
            idbasin(ii) = 1;
        end
    end
    if idbasin(ii)==0 %no vertex was found in basin 1
        for vertex=1:3
            if any(iid2==md.mesh.elements(ii,vertex)) %one vertex in basin 2
                idbasin(ii) = 2;
            end
        end
    end
    if idbasin(ii)==0 %no vertex was found in basin 1 and 2
        idbasin(ii) = 3;
    end
end

%SMB parameters
numparams                  = 2;
numbreaks                  = 1;
intercept                  = [0.5,1.0;1.0,0.6;,2.0,3.0]; %intercept values of SMB in basins [m ice eq./yr]
trendlin                   = [0.0,0.0;0.01,0.001;-0.01,0]; %trend values of SMB in basins [m ice eq./yr^2]
polynomialparams           = cat(3,intercept,trendlin);
datebreaks                 = [3;3;3];

md.timestepping.start_time = 0;
md.timestepping.time_step  = 1/12;
md.timestepping.final_time = 2;
md.smb                     = SMBarma();
md.smb.num_basins          = 3; %number of basins
md.smb.basin_id            = idbasin; %prescribe basin ID number to elements
md.smb.num_params          = numparams; %number of parameters in the polynomial
md.smb.num_breaks          = numbreaks; %number of breakpoints
md.smb.polynomialparams    = polynomialparams;
md.smb.datebreaks          = datebreaks;
md.smb.ar_order            = 4;
md.smb.ma_order            = 1;
md.smb.arma_timestep       = 2.0; %timestep of the ARMA model [yr]
md.smb.arlag_coefs         = [[0.2,0.1,0.05,0.01];[0.4,0.2,-0.2,0.1];[0.4,-0.4,0.1,-0.1]];
md.smb.malag_coefs         = [1.0;0;0.2];
lm0                        = [1e-4*[1,-0.1,-1];1e-6*[1,-0.1,-1];1e-5*[1,-0.1,-1]];
lm1                        = [1e-4*[2,-0.2,-2];1e-6*[2,-0.2,-2];1e-5*[2,-0.2,-2]];
lm2                        = [1e-4*[3,-0.3,-3];1e-6*[3,-0.3,-3];1e-5*[3,-0.3,-3]];
lm3                        = [1e-4*[4,-0.4,-4];1e-6*[4,-0.4,-4];1e-5*[4,-0.4,-4]];
lm4                        = [1e-4*[5,-0.5,-5];1e-6*[5,-0.5,-5];1e-5*[5,-0.5,-5]];
lm5                        = [1e-4*[6,-0.6,-6];1e-6*[6,-0.6,-6];1e-5*[6,-0.6,-6]];
lm6                        = [1e-4*[7,-0.7,-7];1e-6*[7,-0.7,-7];1e-5*[7,-0.7,-7]];
lm7                        = [1e-4*[8,-0.8,-8];1e-6*[8,-0.8,-8];1e-5*[8,-0.8,-8]];
lm8                        = [1e-4*[9,-0.9,-9];1e-6*[9,-0.9,-9];1e-5*[9,-0.9,-9]];
lm9                        = [1e-4*[10,-1,-10];1e-6*[10,-1.0,-10];1e-5*[10,-1.0,-10]];
lm10                       = [1e-4*[11,-1.1,-11];1e-6*[11,-1.1,-11];1e-5*[11,-1.1,-11]];
lm11                       = [1e-4*[12,-1.2,-12];1e-6*[12,-1.2,-12];1e-5*[12,-1.2,-12]];
md.smb.lapserates          = cat(3,lm0,lm1,lm2,lm3,lm4,lm5,lm6,lm7,lm8,lm9,lm10,lm11);
md.smb.elevationbins       = repmat([100,300;200,400;250,450],1,1,12);

%Stochastic forcing
md.stochasticforcing.isstochasticforcing = 1;
md.stochasticforcing.fields              = [{'SMBarma'}];
md.stochasticforcing.covariance          = [[0.15 0.08 -0.02];[0.08 0.12 -0.05];[-0.02 -0.05 0.1]]; %global covariance among- and between-fields
md.stochasticforcing.randomflag          = 0; %fixed random seeds
md.stochasticforcing.stochastictimestep  = 1.0;

md=solve(md,'Transient');

%Fields and tolerances to track changes
field_names     ={'Vx1','Vy1','Vel1','Thickness1','Volume1','SmbMassBalance1','Vx2','Vy2','Vel2','Thickness2','Volume2','SmbMassBalance2','Vx3','Vy3','Vel3','Thickness3','Volume3','SmbMassBalance3'};
field_tolerances={1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,...
						1e-13,1e-13,1e-13,1e-13,1e-13,1e-13...
						1e-12,1e-12,1e-12,1e-12,1e-12,1e-12};
field_values={...
	(md.results.TransientSolution(1).Vx),...
	(md.results.TransientSolution(1).Vy),...
	(md.results.TransientSolution(1).Vel),...
	(md.results.TransientSolution(1).Thickness),...
	(md.results.TransientSolution(1).IceVolume),...
	(md.results.TransientSolution(1).SmbMassBalance),...
	(md.results.TransientSolution(12).Vx),...
	(md.results.TransientSolution(12).Vy),...
	(md.results.TransientSolution(12).Vel),...
	(md.results.TransientSolution(12).Thickness),...
	(md.results.TransientSolution(12).IceVolume),...
	(md.results.TransientSolution(12).SmbMassBalance),...
	(md.results.TransientSolution(24).Vx),...
	(md.results.TransientSolution(24).Vy),...
	(md.results.TransientSolution(24).Vel),...
	(md.results.TransientSolution(24).Thickness),...
	(md.results.TransientSolution(24).IceVolume),...
	(md.results.TransientSolution(24).SmbMassBalance),...
	};
