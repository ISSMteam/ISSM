function md=SetMOLHOBC(md)
%SETMOLHOBC - Create the boundary conditions for stressbalance for MOLHO: VxBase, VyBase, VxShear, VyShear
%
%   Usage:
%      md=SetMOLHOBC(md)
%


%node on Dirichlet
if md.flowequation.isMOLHO
	md.stressbalance.spcvx_base=md.stressbalance.spcvx;
	md.stressbalance.spcvy_base=md.stressbalance.spcvy;

	md.stressbalance.spcvx_shear=NaN*ones(size(md.stressbalance.spcvx_base));
	md.stressbalance.spcvy_shear=NaN*ones(size(md.stressbalance.spcvy_base));
end
