%Test Name: SquareShelfConstrainedTherTranAdolc
md=triangle(model(),'../Exp/Square.exp',180000.);
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareShelfConstrained.par');
md=extrude(md,3,1.);
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',1);
md.transient.isstressbalance=0;
md.transient.ismasstransport=0;
md.transient.issmb=1;
md.transient.isthermal=1;
md.transient.isgroundingline=0;
md.autodiff.isautodiff=true;
md.verbose=verbose('autodiff',true);
md.toolkits.DefaultAnalysis=issmgslsolver();
md=solve(md,'Transient');

%Fields and tolerances to track changes
field_names     ={'Temperature','BasalforcingsGroundediceMeltingRate'};
field_tolerances={1e-13,1e-13};
field_values={...
	(md.results.TransientSolution(1).Temperature),...
	(md.results.TransientSolution(1).BasalforcingsGroundediceMeltingRate),...
	};
