!vim runme.m
runme('id',[101]);
md.mesh.numberofelements
md=solve(md,TransientSolutionEnum);
md=solve(md,StressbalanceSolutionEnum);
plotmodel(md,'data',md.results.StressbalanceSolution.Vel)
plotmodel(md,'data',)
