%Test Name: SquareShelfLevelsetKillberg
md=triangle(model(),'../Exp/Square.exp',50000.);
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareSheetShelf.par');
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',3);

x = md.mesh.x;
xmin = min(x);
xmax = max(x);
Lx = (xmax-xmin);
alpha = 2./3.;
md.mask.ice_levelset = ((x - alpha*Lx)>0) - ((x - alpha*Lx)<0);

% a very special border case that fails the previous version of killberg 20220427
icebergid = 59;
ids = [141,144,139, 58, 263];
md.mask.ice_levelset(ids) = 1;
md.mask.ocean_levelset = -md.mask.ice_levelset;
md.mask.ocean_levelset(142) = 1;
md.mask.ocean_levelset(141) = 1;
md.mask.ocean_levelset(143) = 1;
md.mask.ocean_levelset(144) = 1;
md.mask.ocean_levelset(139) = 1;
md.mask.ice_levelset(icebergid) = -1;

md.levelset.kill_icebergs=1;

md.timestepping.time_step=10;
md.timestepping.final_time=30;

%Transient
md.transient.isstressbalance=0;
md.transient.ismasstransport=1;
md.transient.issmb=1;
md.transient.isthermal=0;
md.transient.isgroundingline=0;
md.transient.ismovingfront=1;

md.calving.calvingrate=zeros(md.mesh.numberofvertices,1);
md.frontalforcings.meltingrate=10000*ones(md.mesh.numberofvertices,1);
md.levelset.spclevelset=md.mask.ice_levelset;

md=solve(md,'Transient');

%Fields and tolerances to track changes
field_names     ={'MaskIceLevelsetwithIceberg','MaskIceLevelset1',...
	'MaskIceLevelset2',...
	'MaskIceLevelset3'};
field_tolerances={1e-11,1e-11,1e-11,1e-11,1e-11,1e-11,1e-11,...
	2e-11,2e-11,2e-11,1e-11,1e-11,1e-11,5e-11,...
	2e-11,2e-11,2e-11,1e-11,1e-11,1e-11,5e-11};
field_values={...
	md.mask.ice_levelset,...
	md.results.TransientSolution(1).MaskIceLevelset,...
	md.results.TransientSolution(2).MaskIceLevelset,...
	md.results.TransientSolution(3).MaskIceLevelset,...
	};

