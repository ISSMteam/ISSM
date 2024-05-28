%Test Name: SquareSheetTherSteaPlume
md=triangle(model(),'../Exp/Square.exp',180000.);
md=setmask(md,'','');
md=parameterize(md,'../Par/SquareSheetConstrained.par');
md.basalforcings=plumebasalforcings(md);
md.basalforcings.floatingice_melting_rate=zeros(md.mesh.numberofvertices,1);
md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
md.basalforcings.plumex=500000;
md.basalforcings.plumey=500000;
md=extrude(md,3,1.);
md=setflowequation(md,'SSA','all');
md.timestepping.time_step=0.;
md.thermal.requested_outputs={'default','BasalforcingsGeothermalflux'};
md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Thermal');

%Fields and tolerances to track changes
field_names     ={'Temperature','BasalforcingsGroundediceMeltingRate','BasalforcingsGeothermalflux'};
field_tolerances={1e-13,1e-8,1e-13};
field_values={...
	(md.results.ThermalSolution.Temperature),...
	(md.results.ThermalSolution.BasalforcingsGroundediceMeltingRate),...
	(md.results.ThermalSolution.BasalforcingsGeothermalflux),...
	};
