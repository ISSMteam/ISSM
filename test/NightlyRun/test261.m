%Test Name: SquareShelfConstrainedTranEnhanced
md=triangle(model(),'../Exp/Square.exp',180000.);
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareShelfConstrained.par');
md=extrude(md,3,1.);
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',3);
md.materials=matenhancedice();
md.materials.rheology_B=3.15e8*ones(md.mesh.numberofvertices,1);
md.materials.rheology_n=3*ones(md.mesh.numberofelements,1);
md.materials.rheology_E=ones(md.mesh.numberofvertices,1);
md.transient.isstressbalance=1;
md.transient.ismasstransport=0;
md.transient.issmb=1;
md.transient.isthermal=1;
md.transient.isgroundingline=0;
md=solve(md,'Transient');

%Fields and tolerances to track changes
field_names     ={'Vx','Vy','Vel','Temperature','BasalforcingsGroundediceMeltingRate'};
field_tolerances={1e-13,1e-13,1e-13,1e-13,1e-13};
field_values={...
	(md.results.TransientSolution(1).Vx),...
	(md.results.TransientSolution(1).Vy),...
	(md.results.TransientSolution(1).Vel),...
	(md.results.TransientSolution(1).Temperature),...
	(md.results.TransientSolution(1).BasalforcingsGroundediceMeltingRate),...
	};
