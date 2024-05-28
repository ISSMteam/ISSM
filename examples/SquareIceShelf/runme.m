md=model;
md=triangle(md,'DomainOutline.exp',100000);
md=setmask(md,'all','');
md=parameterize(md,'Square.par');
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname,'np',2);
md=solve(md,'Stressbalance');
