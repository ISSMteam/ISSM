%Test Name: 79NorthBedSlop2d
md=triangle(model(),'../Exp/79North.exp',10000.);
md=setmask(md,'../Exp/79NorthShelf.exp','');
md=parameterize(md,'../Par/79North.par');
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
