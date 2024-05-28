//Test Name: SquareShelfConstrainedSurfSlop2d
var md = new model();
triangle(md,square[0],150000.);
setmask(md,'all','');
parameterize(md);
setflowequation(md,'SSA','all');
//md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'SurfaceSlope');

//Fields and tolerances to track changes
field_names     =['SurfaceSlopeX','SurfaceSlopeY'];
field_tolerances=[1e-13,1e-13];
field_values=[
	(md.results.SurfaceSlopeSolution[0].SurfaceSlopeX),
	(md.results.SurfaceSlopeSolution[0].SurfaceSlopeY),
	];
