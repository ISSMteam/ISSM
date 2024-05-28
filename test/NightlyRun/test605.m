%Test Name: 79NorthSurfSlop3d
md=triangle(model(),'../Exp/79North.exp',10000.);
md=setmask(md,'../Exp/79NorthShelf.exp','');
md=parameterize(md,'../Par/79North.par');
md=extrude(md,5,1.5);
md=setflowequation(md,'HO','all');
md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'SurfaceSlope');

%Fields and tolerances to track changes
field_names     ={'SurfaceSlopeX','SurfaceSlopeY'};
field_tolerances={1e-13,1e-13};
field_values={...
	(md.results.SurfaceSlopeSolution.SurfaceSlopeX),...
	(md.results.SurfaceSlopeSolution.SurfaceSlopeY),...
	};
