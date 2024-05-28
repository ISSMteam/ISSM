%Test Name: SquareShelfDamageEvolutionSSA2dPralong
md=triangle(model(),'../Exp/Square.exp',150000.);
md=setmask(md,'all','');
md.materials=matdamageice();
md=parameterize(md,'../Par/SquareShelf.par');
md.damage.isdamage=1;
md.damage.D=0.1*ones(md.mesh.numberofvertices,1);
md.damage.spcdamage=NaN*ones(md.mesh.numberofvertices,1);
md.damage.law=1;

md.damage.c1=1.e-11;
md.damage.c2=0.4;
md.damage.c3=1.e-3;
md.damage.healing=0.4;
md.damage.stress_threshold=1e5;
md.damage.stabilization=1;

md.damage.requested_outputs={'default','DamageF'};

md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'DamageEvolution');

field_names={'D','F'};
field_tolerances={1e-13,1e-13};
field_values={...
		(md.results.DamageEvolutionSolution.DamageDbar),...
		(md.results.DamageEvolutionSolution.DamageF),...
	};
