%Test Name: SquareShelfConstrainedEnthalpyStea
md=triangle(model(),'../Exp/Square.exp',180000.);
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareShelfConstrained.par');
md=extrude(md,3,1.);
md=setflowequation(md,'SSA','all');
md.timestepping.time_step=0;
md.initialization.waterfraction=zeros(md.mesh.numberofvertices,1);
md.initialization.watercolumn=zeros(md.mesh.numberofvertices,1);
md.thermal.isenthalpy = 1;
md.thermal.isdynamicbasalspc = 1;

md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Thermal');

%Fields and tolerances to track changes
field_names     ={'Enthalpy','Waterfraction','Temperature'};
field_tolerances={1e-13,3e-10,1e-13};
field_values={...
	(md.results.ThermalSolution.Enthalpy),...
	(md.results.ThermalSolution.Waterfraction),...
	(md.results.ThermalSolution.Temperature),...
	};
