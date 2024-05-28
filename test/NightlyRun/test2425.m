%Test Name: SquareSheetShelfGroundingLine2dSoft
md=triangle(model(),'../Exp/Square.exp',150000.);
md=setmask(md,'../Exp/SquareShelf.exp','');
md=parameterize(md,'../Par/SquareSheetShelf.par');
md=setflowequation(md,'SSA','all');
md.initialization.vx(:)=0.;
md.initialization.vy(:)=0.;
md.geometry.base=-700.-(md.mesh.y-500000.)/1000.;
md.geometry.bed =-700.-(md.mesh.y-500000.)/1000.;
md.geometry.thickness(:)=1300.;
md.geometry.surface=md.geometry.base+md.geometry.thickness;
md.transient.isstressbalance=1;
md.transient.isgroundingline=1;
md.groundingline.migration='AggressiveMigration';

md.timestepping.time_step=.1;
md.timestepping.final_time=1;

md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Transient');
vel1=md.results.TransientSolution(end).Vel;


%get same results with offset in bed and sea level: 
md.geometry.base=-700.-(md.mesh.y-500000.)/1000.;
md.geometry.bed =-700.-(md.mesh.y-500000.)/1000.;
md.geometry.thickness(:)=1300.;
md.geometry.surface=md.geometry.base+md.geometry.thickness;

md.geometry.base=md.geometry.base+1000;
md.geometry.bed=md.geometry.bed+1000;
md.geometry.surface=md.geometry.surface+1000;
md.solidearth.initialsealevel=1000*ones(md.mesh.numberofvertices,1);

md=solve(md,'Transient','checkconsistency','no');
vel2=md.results.TransientSolution(end).Vel;

%Fields and tolerances to track changes
field_names     ={'Vel','Veloffset'};
field_tolerances={1e-13,1e-13};
field_values={vel1,vel2};
