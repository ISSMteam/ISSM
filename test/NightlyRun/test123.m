%Test Name: SquareShelfConstrainedTranMisfitSurface
md=triangle(model(),'../Exp/Square.exp',180000.);
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareShelfConstrained.par');
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',3);
%md.debug.valgrind=1;

fake_surface=[[md.geometry.surface+100;1.1],...
[md.geometry.surface+200;2.1],...
[md.geometry.surface+300;2.5]];

md.transient.requested_outputs={'default','SurfaceMisfit'};
md.outputdefinition.definitions={misfit('name','SurfaceMisfit', 'definitionstring','Outputdefinition1','model_string','Surface','observation_string','SurfaceObservation','observation',fake_surface,'timeinterpolation','nearestneighbor','weights',ones(md.mesh.numberofvertices,1),'weights_string','WeightsSurfaceObservation')};

md=solve(md,'Transient');

%Fields and tolerances to track changes
field_names     ={'SurfaceMisfitFirstStep','SurfaceMisfitSecondStep','SurfaceMisfitThirdStep'};
field_tolerances={1e-13,1e-13,1e-13};
field_values={...
	(md.results.TransientSolution(1).SurfaceMisfit),...
	(md.results.TransientSolution(2).SurfaceMisfit),...
	(md.results.TransientSolution(3).SurfaceMisfit)
	};
