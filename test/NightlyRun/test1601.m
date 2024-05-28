%Test Name: SquareShelfSSA2dRotation
md=triangle(model(),'../Exp/Square.exp',150000.);
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareShelf.par');
md=setflowequation(md,'SSA','all');
md.stressbalance.spcvx(find(md.mesh.y>0.))=NaN;
md.initialization.vx(:)=0.;
md.initialization.vy(:)=0.;
md.initialization.vel(:)=0.;

md.cluster=generic('name',oshostname(),'np',2);
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

plotmodel(md,'data',vel0,'data',vel1,'data',vel1-vel0,'title','Cartesian CS','title','Rotated CS','title','difference')
disp(['Error between Cartesian and rotated CS: ' num2str(max(abs(vel0-vel1))/(max(abs(vel0))+eps)) ]);

%Now, put CS back to normal except on the side where the spc are applied
pos=find(x==0. | x==1000000.);
md.stressbalance.referential(:)=NaN;
md.stressbalance.referential(pos,1:3)=repmat([cos(theta),sin(theta),0],size(pos,1),1);
md.stressbalance.referential(pos,4:6)=repmat([0,0,1],size(pos,1),1);
md=solve(md,'Stressbalance');
vel2=md.results.StressbalanceSolution.Vel;

plotmodel(md,'data',vel0,'data',vel2,'data',vel2-vel0,'title','Cartesian CS','title','Rotated CS','title','difference')
disp(['Error between Cartesian and rotated CS: ' num2str(max(abs(vel0-vel2))/(max(abs(vel0))+eps)) ]);

%Fields and tolerances to track changes
field_names     ={'vel1','vel2'};
field_tolerances={1e-11,1e-11};
field_values={...
	vel1, ...
	vel2, ...
	};
