//Test Name: SquareShelfConstrainedBedSlop3d
var md = new model();
triangle(md,square[0],150000.);
setmask(md,'all','');
parameterize(md);
md.extrude(md,5,1.);
setflowequation(md,'SSA','all');
//md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'BedSlope');

//Fields and tolerances to track changes
field_names     =['BedSlopeX','BedSlopeY'];
field_tolerances=[1e-13,1e-13];
field_values=[
	(md.results.BedSlopeSolution[0].BedSlopeX),
	(md.results.BedSlopeSolution[0].BedSlopeY),
	];
