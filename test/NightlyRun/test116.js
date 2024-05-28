//Test Name: SquareShelfConstrainedBalThic2d
var md = new model();
triangle(md,square[0],150000.);
setmask(md,'all','');
parameterize(md);
//Add boundary conditions on thickness on the border
pos=[];
for (var i = 0; i < md.mesh.vertexonboundary.length; ++i) {
    if ((md.mesh.vertexonboundary[i] !== 0)) {
            pos.push(i);
    };
}

for (var i = 0; i < pos.length; ++i) {
    md.balancethickness.spcthickness[pos[i]] = md.geometry.thickness[pos[i]];
}
setflowequation(md,'SSA','all');
//md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Balancethickness');

//Fields and tolerances to track changes
field_names     =['Thickness'];
field_tolerances=[1e-13];
field_values=[
	(md.results.BalancethicknessSolution[0].Thickness),
	];
