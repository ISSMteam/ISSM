%Test Name: EISMINTTran2
%Test on the stressbalance model and the masstransport in 2d
printingflag=false;

%tests 3 and 4: using Glen's flow law
md=model();
md=triangle(md,'../Exp/SquareEISMINT.exp',3550.);
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareEISMINT.par');
md=setflowequation(md,'SSA','all'); %SSA's model and 2d

%Impose a non zero velocity on the upper boundary condition (y=max(y))
pos=find(md.mesh.y==max(md.mesh.y));
md.stressbalance.spcvy(pos)=400.*(((md.mesh.x(pos)-100000.)/25000.).^2-ones(size(pos,1),1)).*heaviside((1.+eps)*ones(size(pos,1),1)-((md.mesh.x(pos)-100000.)/25000.).^2);

%Compute solution for SSA's model 
md.cluster=generic('name',oshostname(),'np',8);
md=solve(md,'Stressbalance');

%plot results
md.initialization.vx=(md.results.StressbalanceSolution.Vx);
md.initialization.vy=(md.results.StressbalanceSolution.Vy);

md.timestepping.time_step=1.;
md.timestepping.final_time=5000.;
md.masstransport.stabilization=1;
md=solve(md,'Transient');

plotmodel(md,'data',(md.results.TransientSolution(end).Vx))
if printingflag,
	set(gcf,'Color','w')
	printmodel('eisminttrans2vx','png','margin','on','marginsize',25,'frame','off','resolution',2,'hardcopy','off');
	system(['mv eisminttrans2vx.png ' ISSM_DIR '/website/doc_pdf/validation/Images/EISMINT/IceShelf ']);
end

plotmodel(md,'data',(md.results.TransientSolution(end).Vy))
if printingflag,
	set(gcf,'Color','w')
	printmodel('eisminttrans2vy','png','margin','on','marginsize',25,'frame','off','resolution',2,'hardcopy','off');
	system(['mv eisminttrans2vy.png ' ISSM_DIR '/website/doc_pdf/validation/Images/EISMINT/IceShelf ']);
end

plotmodel(md,'data',(md.results.TransientSolution(end).Thickness))
if printingflag,
	set(gcf,'Color','w')
	printmodel('eisminttrans2thickness','png','margin','on','marginsize',25,'frame','off','resolution',2,'hardcopy','off');
	system(['mv eisminttrans2thickness.png ' ISSM_DIR '/website/doc_pdf/validation/Images/EISMINT/IceShelf ']);
end

%Fields and tolerances to track changes
field_names     ={'Vx','Vy','Thickness'};
field_tolerances={1e-13,1e-13,1e-13};
field_values={...
	(md.results.TransientSolution(end).Vx), ...
	(md.results.TransientSolution(end).Vy), ...
	(md.results.TransientSolution(end).Thickness), ...
	};
