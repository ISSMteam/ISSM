%Test Name: SquareSheetShelfStressHOHigherOrder
md=triangle(model(),'../Exp/Square.exp',180000.);
md=setmask(md,'../Exp/SquareShelf.exp','');
md=parameterize(md,'../Par/SquareSheetShelf.par');
md=extrude(md,5,1.);
md=setflowequation(md,'HO','all');
md.cluster=generic('name',oshostname(),'np',3);

field_names={};
field_tolerances={};
field_values={};
for i={'P1bubble','P1bubblecondensed','P1xP2','P2xP1','P2','P1xP3','P2xP4'}
	md.flowequation.fe_HO=i{1};
	md=solve(md,'Stressbalance');
	field_names     ={field_names{:},['Vx' i{1}],['Vy' i{1}],['Vz' i{1}],['Vel' i{1}],['Pressure' i{1}]};
	field_tolerances={field_tolerances{:},1e-07,6e-08,6e-08,6e-08,3e-13};
	field_values={field_values{:},...
	(md.results.StressbalanceSolution.Vx),...
	(md.results.StressbalanceSolution.Vy),...
	(md.results.StressbalanceSolution.Vz),...
	(md.results.StressbalanceSolution.Vel),...
	(md.results.StressbalanceSolution.Pressure),...
	};
end
