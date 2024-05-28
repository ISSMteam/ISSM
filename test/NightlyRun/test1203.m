%Test Name: EISMINTStress2
%Test on the stressbalance model and the masstransport in 2d
printingflag=false;

%test 5 and 6: 
md=model();
md=triangle(md,'../Exp/SquareEISMINT.exp',5100.); %test3
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareEISMINT.par');
md=setflowequation(md,'SSA','all'); %SSA's model and 2d

%Impose a non zero velocity on the upper boundary condition (y=max(y))
pos=find(md.mesh.y==max(md.mesh.y));
md.stressbalance.spcvy(pos)=400.*(((md.mesh.x(pos)-100000.)/25000.).^2-ones(size(pos,1),1)).*heaviside((1.+eps)*ones(size(pos,1),1)-((md.mesh.x(pos)-100000.)/25000.).^2);

%Compute solution for SSA's model 
md.cluster=generic('name',oshostname(),'np',8);
md=solve(md,'Stressbalance');

vx=(md.results.StressbalanceSolution.Vx);
vy=(md.results.StressbalanceSolution.Vy);

%plot results
plotmodel(md,'data',vx,'contourlevels',{0,20,40,60,80,100,-20,-40,-60,-80,-100},...
	'contourcolor','k')
if printingflag,
	set(gcf,'Color','w')
	printmodel('eismintdiag2vx','png','margin','on','marginsize',25,'frame','off','resolution',2,'hardcopy','off');
	system(['mv eismintdiag2vx.png ' ISSM_DIR '/website/doc_pdf/validation/Images/EISMINT/IceShelf ']);
end
plotmodel(md,'data',vy,'contourlevels',{-100,-200,-300,-400,-500,-600,-700,-800,-900,-1000},...
	'contourcolor','k')
if printingflag,
	set(gcf,'Color','w')
	printmodel('eismintdiag2vy','png','margin','on','marginsize',25,'frame','off','resolution',2,'hardcopy','off');
	system(['mv eismintdiag2vy.png ' ISSM_DIR '/website/doc_pdf/validation/Images/EISMINT/IceShelf ']);
end

%Fields and tolerances to track changes
field_names     ={'Vx','Vy'};
field_tolerances={1e-13,1e-13};
field_values={...
	vx, ...
	vy, ...
	};
