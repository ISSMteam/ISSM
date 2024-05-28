//Test Name: SquareShelfConstrained
var md = new model();
function zeros(...args) {
	var array = [];
	for (var i = 0; i < args[0]; ++i) {
		array.push(args.length == 1 ? 0 : zeros(args.slice(1)));
	}
	return array;
}
var md = new model();
triangle(md,square[0],150000.);
setmask(md,'all','');
parameterize(md);
setflowequation(md,'SSA','all');
//md.cluster=generic('name',oshostname(),'np',3);

//redo the parameter file for this special shelf. 
//constant thickness, constrained (vy=0) flow into an icefront, 
//from 0 m/yr at the grounding line.

//tighten
md.stressbalance.restol=Math.pow(10,-4);

//needed later
ymin=Math.min.apply(null, md.mesh.y);
ymax=Math.max.apply(null, md.mesh.y);
xmin=Math.min.apply(null, md.mesh.x);
xmax=Math.max.apply(null, md.mesh.x);

di=md.materials.rho_ice/md.materials.rho_water;

h=1000;
md.geometry.thickness=h*ones(md.mesh.numberofvertices,1);
md.geometry.base=-md.materials.rho_ice/md.materials.rho_water*md.geometry.thickness;
md.geometry.surface=md.geometry.base+md.geometry.thickness;

//Initial velocity and pressure
md.initialization.vx=zeros(md.mesh.numberofvertices,1);
md.initialization.vy=zeros(md.mesh.numberofvertices,1);
md.initialization.vz=zeros(md.mesh.numberofvertices,1);
md.initialization.pressure=zeros(md.mesh.numberofvertices,1);

//Materials
md.initialization.temperature=(273-20)*ones(md.mesh.numberofvertices,1);
md.materials.rheology_B=paterson(md.initialization.temperature);
md.materials.rheology_n=3*ones(md.mesh.numberofelements,1);

//Boundary conditions:
md.stressbalance.spcvx=NaN*ones(md.mesh.numberofvertices,1);
md.stressbalance.spcvy=NaN*ones(md.mesh.numberofvertices,1);
md.stressbalance.spcvz=NaN*ones(md.mesh.numberofvertices,1);

//constrain flanks to 0 normal velocity
pos=[];
for (var i = 0; i < md.mesh.x.length; ++i) {
    if (md.mesh.x[i]==xmin || md.mesh.x[i]==xmax) {
            pos.push(i);
    };
}

for (var i = 0; i < pos.length; ++i) {
    md.stressbalance.spcvx[pos[i]] = 0;
    md.stressbalance.spcvz[pos[i]] = NaN;
}

//constrain grounding line to 0 velocity
pos=[];
for (var i = 0; i < md.mesh.y.length; ++i) {
    if (md.mesh.y[i]==ymin) {
            pos.push(i);
    };
}
for (var i = 0; i < pos.length; ++i) {
    md.stressbalance.spcvx[pos[i]] = 0;
    md.stressbalance.spcvy[pos[i]] = 0;
}

//icefront
nodeonicefront=zeros(md.mesh.numberofvertices,1);
for (var i = 0; i < md.mesh.y.length; ++i) {
    if (md.mesh.y[i]==ymax) {
        nodeonicefront[i] = 1;
    }
}
md.mask.ice_levelset = nodeonicefront.map(function(x) { return -1 + x; });

md=solve(md,'Stressbalance');

//create analytical solution: strain rate is constant = ((rho_ice*g*h)/4B)Math.pow(,3) (Paterson, 4th Edition, page 292.
//ey_c=(md.materials.rho_ice*md.constants.g*(1-di)*md.geometry.thickness./(4*md.materials.rheology_B)).Math.pow(,3);
//vy_c=ey_c.*md.mesh.y*md.constants.yts;

//Fields and tolerances to track changes
field_names     =['Vy'];
field_tolerances=[1e-13];
field_values=[(md.results.StressbalanceSolution[0].Vy)];
