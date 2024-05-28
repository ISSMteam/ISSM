%Test Name: SquareShelfStressHOHigherOrder
md=triangle(model(),'../Exp/Square.exp',150000.);
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareShelf.par');
md=extrude(md,3,2.);
md=setflowequation(md,'HO','all');
md.cluster=generic('name',oshostname(),'np',3);

field_names={};
field_tolerances={};
field_values={};
for i={'P1bubble','P1bubblecondensed','P1xP2','P2xP1','P2','P1xP3','P2xP4'}
	md.flowequation.fe_HO=i{1};
	md=solve(md,'Stressbalance');
	field_names     ={field_names{:},['Vx' i{1}],['Vy' i{1}],['Vz' i{1}],['Vel' i{1}],['Pressure' i{1}]};
	field_tolerances={field_tolerances{:},6.7e-08,5e-08,2e-08,5e-08,1e-13};
	field_values={field_values{:},...
	(md.results.StressbalanceSolution.Vx),...
	(md.results.StressbalanceSolution.Vy),...
	(md.results.StressbalanceSolution.Vz),...
	(md.results.StressbalanceSolution.Vel),...
	(md.results.StressbalanceSolution.Pressure),...
	};
end
