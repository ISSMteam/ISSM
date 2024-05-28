//Test Name: SquareShelfConstrainedEnthalpyTran
function zeros(...args) {
	var array = [];
	for (var i = 0; i < args[0]; ++i) {
		array.push(args.length == 1 ? 0 : zeros(args.slice(1)));
	}
	return array;
}
var md = new model();
triangle(md,square[0],180000.);
setmask(md,'all','');
parameterize(md);
md.extrude(md,3,1.);
setflowequation(md,'SSA','all');
//md.cluster=generic('name',oshostname(),'np',3);
md.initialization.waterfraction=zeros(md.mesh.numberofvertices,1);
md.initialization.watercolumn=zeros(md.mesh.numberofvertices,1);
md.transient.isstressbalance=0;
md.transient.ismasstransport=0;
md.transient.issmb=1;
md.transient.isthermal=1;
md.transient.isgroundingline=0;
md.thermal.isenthalpy=1;
md.thermal.isdynamicbasalspc=1;
md=solve(md,'Transient');

//Fields and tolerances to track changes
field_names     =['Enthalpy1','Waterfraction1','Temperature1',
	'Enthalpy2','Waterfraction2','Temperature2',
	'Enthalpy3','Waterfraction3','Temperature3'];
field_tolerances=[1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-9,1e-13];
field_values=[
	(md.results.TransientSolution[0](1).Enthalpy),
	(md.results.TransientSolution[0](1).Waterfraction),
	(md.results.TransientSolution[0](1).Temperature),
	(md.results.TransientSolution[0](2).Enthalpy),
	(md.results.TransientSolution[0](2).Waterfraction),
	(md.results.TransientSolution[0](2).Temperature),
	(md.results.TransientSolution[0](3).Enthalpy),
	(md.results.TransientSolution[0](3).Waterfraction),
	(md.results.TransientSolution[0](3).Temperature),
	];
