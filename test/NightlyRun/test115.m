%Test Name: SquareShelfConstrainedBedSlop3d
md=triangle(model(),'../Exp/Square.exp',150000.);
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareShelfConstrained.par');
md=extrude(md,5,1.);
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'BedSlope');

%Fields and tolerances to track changes
field_names     ={'BedSlopeX','BedSlopeY'};
field_tolerances={1e-13,1e-13};
field_values={...
	(md.results.BedSlopeSolution.BedSlopeX),...
	(md.results.BedSlopeSolution.BedSlopeY),...
	};
