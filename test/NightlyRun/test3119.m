%Test Name: ReverseScalarDriver
%reverse scalar driver in ADOLC, using the test3009 setup, equivalent to test109 setup.
md=triangle(model(),'../Exp/Square.exp',100000.);
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareShelfConstrained.par');
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',1);
md.toolkits.DefaultAnalysis=issmgslsolver();

md.autodiff.isautodiff=true;
md.verbose.autodiff=true;

%first run scalar reverse mode: 
md.autodiff.independents={independent('name','md.geometry.thickness','type','vertex','nods',md.mesh.numberofvertices)};
md.autodiff.dependents={dependent('name','MaxVel','type','scalar','fos_reverse_index',1)};
md.autodiff.driver='fos_reverse';

md=solve(md,'Transient');

%recover jacobian: 
jac_reverse=md.results.TransientSolution(1).AutodiffJacobian;

%Fields and tolerances to track changes
field_names     ={'Jac Reverse'};
field_tolerances={1e-8};
field_values={jac_reverse,};
