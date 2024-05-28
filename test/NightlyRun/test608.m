%Test Name: 79NorthBalThic2d
md=triangle(model(),'../Exp/79North.exp',10000.);
md=setmask(md,'../Exp/79NorthShelf.exp','');
md=parameterize(md,'../Par/79North.par');
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Balancethickness');

%Fields and tolerances to track changes
field_names     ={'Thickness'};
field_tolerances={1e-12};
field_values={...
	(md.results.BalancethicknessSolution.Thickness),...
	};
