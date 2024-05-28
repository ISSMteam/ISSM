//Test Name: SquareShelfConstrainedTranMisfitSurface
var md = new model();
triangle(md,square[0],180000.);
setmask(md,'all','');
parameterize(md);
setflowequation(md,'SSA','all');
//md.cluster=generic('name',oshostname(),'np',3);

fake_surface=[[md.geometry.surface+100,1.1],[md.geometry.surface+200,2.1],[md.geometry.surface+300,2.5]];

md.transient.requested_outputs=['default','SurfaceMisfit'];
md.outputdefinition.definitions=[misfit('name','SurfaceMisfit', 'definitionenum',Outputdefinition1Enum,'model_enum',SurfaceEnum,'observation_enum',SurfaceObservationEnum,'observation',fake_surface,'timeinterpolation','nearestneighbor','weights',ones(md.mesh.numberofvertices,1),'weights_enum',WeightsSurfaceObservationEnum)];

md=solve(md,'Transient');

//Fields and tolerances to track changes
field_names     =['SurfaceMisfitFirstStep','SurfaceMisfitSecondStep','SurfaceMisfitThirdStep'];
field_tolerances=[1e-13,1e-13,1e-13];
field_values=[
	(md.results.TransientSolution[0](1).SurfaceMisfit),
	(md.results.TransientSolution[0](2).SurfaceMisfit),
	(md.results.TransientSolution[0](3).SurfaceMisfit)
	];
