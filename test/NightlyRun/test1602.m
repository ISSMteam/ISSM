%Test Name: SquareSheetShelfHORotation
md=triangle(model(),'../Exp/Square.exp',150000.);
md=setmask(md,'../Exp/SquareShelf.exp','');
md=parameterize(md,'../Par/SquareSheetShelf.par');
md=extrude(md,5,1.);
md=setflowequation(md,'HO','all');
md.stressbalance.spcvx(find(md.mesh.y>0.))=NaN;
md.initialization.vx(:)=0.;
md.initialization.vy(:)=0.;
md.initialization.vel(:)=0.;

md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Stressbalance');
vel0=md.results.StressbalanceSolution.Vel;

theta=30.*pi/180.;
x=md.mesh.x;
y=md.mesh.y;
md.mesh.x=cos(theta)*x-sin(theta)*y;
md.mesh.y=sin(theta)*x+cos(theta)*y;

md.stressbalance.referential(:,1:3)=repmat([cos(theta),sin(theta),0],md.mesh.numberofvertices,1);
md.stressbalance.referential(:,4:6)=repmat([0,0,1],md.mesh.numberofvertices,1);
md=solve(md,'Stressbalance');
vel1=md.results.StressbalanceSolution.Vel;

plotmodel(md,'data',vel0,'data',vel1,'data',vel1-vel0,'title','Cartesian CS','title','Rotated CS','title','difference','view#all',2)
disp(['Error between Cartesian and rotated CS: ' num2str(max(abs(vel0-vel1))/(max(abs(vel0))+eps)) ]);

%Fields and tolerances to track changes
field_names     ={'vel1'};
field_tolerances={1e-9};
field_values={...
	vel1, ...
	};
