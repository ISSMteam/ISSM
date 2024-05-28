%Test Name: SquareSheetShelfSteaEnthalpyRheologiesHO
md=triangle(model(),'../Exp/Square.exp',150000.);
md=setmask(md,'../Exp/SquareShelf.exp','');
md=parameterize(md,'../Par/SquareSheetShelf.par');
md=extrude(md,3,2.);
md=setflowequation(md,'HO','all');
md.cluster=generic('name',oshostname(),'np',3);
md.timestepping.time_step=0.;
md.thermal.isenthalpy=1;
md.initialization.waterfraction=zeros(md.mesh.numberofvertices,1);
md.initialization.watercolumn=zeros(md.mesh.numberofvertices,1);

%Go solve
field_names={};
field_tolerances={};
field_values={};
for i={'LliboutryDuval', 'CuffeyTemperate'}
	disp(' ');
	disp(['====== Testing rheology law: ' i{1} ' =====']);

	md.materials.rheology_law=i{1};
	md=solve(md,'Steadystate');
	field_names     ={field_names{:},['Vx' i{1}],['Vy' i{1}],['Vz' i{1}],['Vel' i{1}],['Pressure' i{1}],...
		['Temperature' i{1}],['Waterfraction' i{1}],['Enthalpy' i{1}]};
	field_tolerances={field_tolerances{:},2e-09,1e-09,1e-09,1e-09,1e-13,2e-10,6e-10,1e-9};
	field_values={field_values{:},...
		(md.results.SteadystateSolution.Vx),...
		(md.results.SteadystateSolution.Vy),...
		(md.results.SteadystateSolution.Vz),...
		(md.results.SteadystateSolution.Vel),...
		(md.results.SteadystateSolution.Pressure),...
		(md.results.SteadystateSolution.Temperature),...
		(md.results.SteadystateSolution.Waterfraction),...
		(md.results.SteadystateSolution.Enthalpy),...
		};
end
