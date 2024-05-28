%Test Name: SquareShelfStressSSA2dEnhanced
md=triangle(model(),'../Exp/Square.exp',150000.);
md=setmask(md,'all','');
md.materials=matenhancedice();
md.materials.rheology_B=3.15e8*ones(md.mesh.numberofvertices,1);
md.materials.rheology_n=3*ones(md.mesh.numberofelements,1);
md.materials.rheology_E=ones(md.mesh.numberofvertices,1);
md=parameterize(md,'../Par/SquareShelf.par');
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Stressbalance');

%Fields and tolerances to track changes
field_names     ={'Vx','Vy','Vel','Pressure'};
field_tolerances={1e-13,1e-13,1e-13,1e-13};
field_values={...
	(md.results.StressbalanceSolution.Vx),...
	(md.results.StressbalanceSolution.Vy),...
	(md.results.StressbalanceSolution.Vel),...
	(md.results.StressbalanceSolution.Pressure),...
	};
