%Test Name: 79NorthMasstransp2dDG
md=triangle(model(),'../Exp/79North.exp',10000.);
md=meshconvert(md);
md=setmask(md,'../Exp/79NorthShelf.exp','');
md=parameterize(md,'../Par/79North.par');
md=setflowequation(md,'SSA','all');
md.masstransport.stabilization=3;
md.masstransport.spcthickness=md.geometry.thickness;
md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Masstransport');

%Fields and tolerances to track changes
field_names     ={'Thickness'};
field_tolerances={1e-13};
field_values={...
	(md.results.MasstransportSolution.Thickness),...
	};
