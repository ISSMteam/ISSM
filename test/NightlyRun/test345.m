%Test Name: SquareSheetConstrainedEnthalpyEffectiveConductivity
md=triangle(model(),'../Exp/Square.exp',180000.);
md=setmask(md,'','');
md=parameterize(md,'../Par/SquareSheetConstrained.par');
md=extrude(md,3,1.);
md=setflowequation(md,'SSA','all');
md.timestepping.time_step=0.;
md.initialization.waterfraction=zeros(md.mesh.numberofvertices,1);
md.initialization.watercolumn=zeros(md.mesh.numberofvertices,1);
md.thermal.isenthalpy = 1;
md.thermal.isdynamicbasalspc = 0;

md.cluster=generic('name',oshostname(),'np',3);
field_names={};
field_tolerances={};
field_values={};
for i=[0 1 2]
	md.materials.effectiveconductivity_averaging=i;
	disp(' ');
	disp(['====== Testing effectiveconductivity_averaging=',num2str(i),' =====']);
	md=solve(md,'Thermal');
	%Fields and tolerances to track changes
	field_names     ={field_names{:},['Enthalpy' i],['Waterfraction' i],['Temperature' i]};
	field_tolerances={field_tolerances{:},1e-13,5e-13,1e-13};
	field_values={field_values{:},...
		(md.results.ThermalSolution.Enthalpy),...
		(md.results.ThermalSolution.Waterfraction),...
		(md.results.ThermalSolution.Temperature),...
		};
end
