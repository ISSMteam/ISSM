//Test Name: SquareShelfConstrainedEnthalpyStea
var md = new model();
function zeros(...args) {
	var array = [];
	for (var i = 0; i < args[0]; ++i) {
		array.push(args.length == 1 ? 0 : zeros(args.slice(1)));
	}
	return array;
}
triangle(md,square[0],180000.);
setmask(md,'all','');
parameterize(md);
md.extrude(md,3,1.);
setflowequation(md,'SSA','all');
md.timestepping.time_step=0;
md.initialization.waterfraction=zeros(md.mesh.numberofvertices,1);
md.initialization.watercolumn=zeros(md.mesh.numberofvertices,1);
md.thermal.isenthalpy = 1;
md.thermal.isdynamicbasalspc = 1;

//md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Thermal');

//Fields and tolerances to track changes
field_names     =['Enthalpy','Waterfraction','Temperature'];
field_tolerances=[1e-13,3e-10,1e-13];
field_values=[
	(md.results.ThermalSolution[0].Enthalpy),
	(md.results.ThermalSolution[0].Waterfraction),
	(md.results.ThermalSolution[0].Temperature),
	];
