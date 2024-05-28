%Test Name: SquareShelfStressSSA2dDamageRift
md=triangle(model(),'../Exp/SquareHole.exp','../Exp/Rifts.exp',50000.);
md=meshprocessrifts(md,'../Exp/Square.exp');
md=setmask(md,'all','');
md.materials=matdamageice();
md=parameterize(md,'../Par/SquareShelf2.par');
md.damage.isdamage=1;
md.damage.D=0.5*ones(md.mesh.numberofvertices,1);
md.damage.spcdamage=NaN*ones(md.mesh.numberofvertices,1);
md=setflowequation(md,'SSA','all');

md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Stressbalance');

%Fields and tolerances to track changes
field_names     ={'Vx','Vy','Vel','Pressure'};
field_tolerances={7e-8,3e-8,3e-8,1e-11};
field_values={...
	   (md.results.StressbalanceSolution.Vx),...
	   (md.results.StressbalanceSolution.Vy),...
	   (md.results.StressbalanceSolution.Vel),...
	   (md.results.StressbalanceSolution.Pressure),...
	   };
