console.log('creating model'); 
var md = new model();

console.log('meshing');
triangle(md,square[0],40000); 

console.log('parameterization');
setmask(md,'all','');
parameterize(md);
setflowequation(md,'SSA','all');
md.verbose.solution=1;  md.verbose.convergence=0;

md=solve(md,'Stressbalance','checkconsistency','no');

console.log(md.results['StressbalanceSolution'][0]['Vel']);
