%Test Name: ISMIPCFS
%This test is a test from the ISMP-HOM Intercomparison project.
%Pattyn and Payne 2006

%L_list={5000.,10000.,20000.,40000.,80000.,160000.};
L_list={80000.};
results={};

for i=1:length(L_list),
	L=L_list{i};  
	md=triangle(model(),['../Exp/Square_' num2str(L) '.exp'],L/10.); %size 3*L 
	md=setmask(md,'',''); %ice sheet test
	md=parameterize(md,'../Par/ISMIPC.par');
	md.friction.coefficient=sqrt(md.constants.yts.*(1000.+1000.*sin(md.mesh.x*2.*pi/L).*sin(md.mesh.y*2.*pi/L)));
	md=extrude(md,10,1.);

	%Add spc on the borders
	pos=find(md.mesh.x==0. | md.mesh.x==max(md.mesh.x) | md.mesh.y==0. | md.mesh.y==max(md.mesh.y));
	md.stressbalance.spcvx(pos)=0.;
	md.stressbalance.spcvy(pos)=0.;
	if(L==5000.),
		md.stressbalance.spcvx(pos)=15.66;
		md.stressbalance.spcvy(pos)=-0.1967;
	elseif(L==10000.),
		md.stressbalance.spcvx(pos)=16.04;
		md.stressbalance.spcvy(pos)=-0.1977;
	elseif(L==20000.),
		md.stressbalance.spcvx(pos)=16.53;
		md.stressbalance.spcvy(pos)=-1.27;
	elseif(L==40000.),
		md.stressbalance.spcvx(pos)=17.23;
		md.stressbalance.spcvy(pos)=-3.17;
	elseif(L==80000.),
		md.stressbalance.spcvx(pos)=16.68;
		md.stressbalance.spcvy(pos)=-2.69;
	elseif(L==160000.),
		md.stressbalance.spcvx(pos)=16.03;
		md.stressbalance.spcvy(pos)=-1.27;
	end

	md=setflowequation(md,'FS','all');

	%Compute the stressbalance
	md.cluster=generic('name',oshostname(),'np',8);
	md=solve(md,'Stressbalance');

	%Plot the results and save them
	vx=(md.results.StressbalanceSolution.Vx);
	vy=(md.results.StressbalanceSolution.Vy);
	vz=(md.results.StressbalanceSolution.Vz);
	results{i}=md.results.StressbalanceSolution;

	plotmodel(md,'data',vx,'data',vy,'data',vz,'layer#all',md.mesh.numberoflayers)
end

%Fields and tolerances to track changes
field_names     ={...
	'Vx80km','Vy80km','Vz80km'
};
field_tolerances={...
	1e-11,2e-12,3e-12,...
};
field_values={};
for i=1:1,
	result=results{i};
	field_values={field_values{:},...
		(result.Vx),...
		(result.Vy),...
		(result.Vz),...
		};
end
