%Test Name: SquareSheetShelfDiadSSA3dDakotaAreaAverage
%test partitioning, and partition averaging
md=triangle(model(),'../Exp/Square.exp',30000.);
md=setmask(md,'../Exp/SquareShelf.exp','');
md=parameterize(md,'../Par/SquareSheetShelf.par');
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',3);

%partitioning
npart=100;
[partition,md]=partitioner(md,'package','chaco','npart',npart);
partition=partition-1;

vector=(1:1:md.mesh.numberofvertices)';
vector_on_partition=AreaAverageOntoPartition(md,vector,partition);
vector_on_nodes=vector_on_partition(partition+1);

field_names     ={'vector_on_nodes'};
field_tolerances={1e-11};
field_values={...
         vector_on_nodes,...
	};
