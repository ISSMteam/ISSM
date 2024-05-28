%Test Name: SquareShelfStressSSA2dDamageUpdate
md=triangle(model(),'../Exp/Square.exp',150000.);
md=setmask(md,'all','');
md.materials=matdamageice();
md=parameterize(md,'../Par/SquareShelf.par');
md.damage.isdamage=1;
md.damage.D=0.*ones(md.mesh.numberofvertices,1);
md.damage.spcdamage=NaN*ones(md.mesh.numberofvertices,1);
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',3);

md.stressbalance.requested_outputs={'default','NewDamage'};
md.damage.stress_threshold=1.3e5;
md.damage.kappa=2.8;

md=solve(md,'Stressbalance');

%Fields and tolerances to track changes
field_names     ={'Vx','Vy','Vel','Pressure','NewDamage'};
field_tolerances={1e-13,1e-13,1e-13,1e-13,1e-13};
field_values={...
	(md.results.StressbalanceSolution.Vx),...
	(md.results.StressbalanceSolution.Vy),...
	(md.results.StressbalanceSolution.Vel),...
	(md.results.StressbalanceSolution.Pressure),...
	(md.results.StressbalanceSolution.NewDamage),...
	};
