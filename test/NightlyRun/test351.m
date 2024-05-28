%Test Name: SquareSheetFrictionTemp
md=triangle(model(),'../Exp/Square.exp',50000.);
md.mesh.x = md.mesh.x/1000;
md.mesh.y = md.mesh.y/1000;
md=setmask(md,'','');
md=parameterize(md,'../Par/SquareSheetConstrained.par');
md=extrude(md,10,3);
md=setflowequation(md,'HO','all');
md.cluster=generic('name',oshostname(),'np',4);

%Use hydroogy coupled friciton law
md.friction=frictiontemp(md.friction);
md.friction.gamma=5;

% Thermal settings
md.thermal.isenthalpy=1;
md.thermal.isdynamicbasalspc=0;
md.thermal.maxiter=5;
md.initialization.waterfraction=zeros(md.mesh.numberofvertices,1);
md.initialization.watercolumn=zeros(md.mesh.numberofvertices,1);

md.basalforcings.geothermalflux=4.2*10^-2*ones(md.mesh.numberofvertices,1);

md.timestepping.time_step=0;
md=solve(md,'Steadystate');

%Fields and tolerances to track changes
field_names     ={'Vx','Vy','Vz'};
field_tolerances={1e-12,1e-12,2e-11};
field_values={...
	(md.results.SteadystateSolution.Vx),...
	(md.results.SteadystateSolution.Vy),...
	(md.results.SteadystateSolution.Vz),...
};
