%Test Name: SquareSheetShelfStressSSAHigherOrder
md=triangle(model(),'../Exp/Square.exp',150000.);
md=setmask(md,'../Exp/SquareShelf.exp','');
md=parameterize(md,'../Par/SquareSheetShelf.par');
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',3);

field_names={};
field_tolerances={};
field_values={};
for i={'P1bubble','P1bubblecondensed','P2'}
	md.flowequation.fe_SSA=i{1};
	md=solve(md,'Stressbalance');
	field_names     ={field_names{:},['Vx' i{1}],['Vy' i{1}],['Vel' i{1}],['Pressure' i{1}]};
	field_tolerances={field_tolerances{:},1e-12,1e-13,1e-13,1e-13};
	field_values={field_values{:},...
		(md.results.StressbalanceSolution.Vx),...
		(md.results.StressbalanceSolution.Vy),...
		(md.results.StressbalanceSolution.Vel),...
		(md.results.StressbalanceSolution.Pressure),...
		};
end
