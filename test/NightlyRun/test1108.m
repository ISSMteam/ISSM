%Test Name: ISMIPDFS
%This test is a test from the ISMP-HOM Intercomparison project.
%Pattyn and Payne 2006

%L_list={5000.,10000.,20000.,40000.,80000.,160000.};
L_list={80000.};
results={};

for i=1:length(L_list),
	L=L_list{i};
	nx=30; %numberof nodes in x direction
	ny=30;
	md=model();
	md=squaremesh(md,L,L,nx,ny);
	md=setmask(md,'',''); %ice sheet test
	md=parameterize(md,'../Par/ISMIPD.par');
	md=extrude(md,10,1.);

	md=setflowequation(md,'HO','all');

	%We need one grid on dirichlet: the 4 corners are set to zero
	md.stressbalance.spcvx=NaN*ones(md.mesh.numberofvertices,1);
	md.stressbalance.spcvy=NaN*ones(md.mesh.numberofvertices,1);
	md.stressbalance.spcvz=NaN*ones(md.mesh.numberofvertices,1);
	
	pos=find(md.mesh.vertexonbase & (md.mesh.x==0. | md.mesh.x==max(md.mesh.x)) & (md.mesh.y==0. | md.mesh.y==max(md.mesh.y)));
	md.stressbalance.spcvx(pos)=0.;
	md.stressbalance.spcvy(pos)=0.;
	md.stressbalance.spcvz(pos)=0.;

	%Create MPCs to have periodic boundary conditions
	posx=find(md.mesh.x==0.);
	posx2=find(md.mesh.x==max(md.mesh.x));

	posy=find(md.mesh.y==0. & md.mesh.x~=0. & md.mesh.x~=max(md.mesh.x)); %Don't take the same nodes two times
	posy2=find(md.mesh.y==max(md.mesh.y) & md.mesh.x~=0. & md.mesh.x~=max(md.mesh.x));

	md.stressbalance.vertex_pairing=[posx,posx2;posy,posy2];

	%Compute the stressbalance
	md.cluster=generic('name',oshostname(),'np',8);
	md=solve(md,'Stressbalance');
	md.stressbalance.reltol=NaN;
	md.stressbalance.abstol=NaN;
	md.stressbalance.vertex_pairing=[];
	%We need one grd on dirichlet: the 4 corners are set to zero
	md.stressbalance.spcvx=NaN*ones(md.mesh.numberofvertices,1);
	md.stressbalance.spcvy=NaN*ones(md.mesh.numberofvertices,1);
	md.stressbalance.spcvz=NaN*ones(md.mesh.numberofvertices,1);
	pos=find(md.mesh.y==0. | md.mesh.x==0. | md.mesh.x==max(md.mesh.x) | md.mesh.y==max(md.mesh.y)); %Don't take the same nodes two times
	md.stressbalance.spcvx(pos)=md.results.StressbalanceSolution.Vx(pos);
	md.stressbalance.spcvy(pos)=md.results.StressbalanceSolution.Vy(pos);
	md=setflowequation(md,'FS','all');
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
	1e-08,1e-07,1e-07,...
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
