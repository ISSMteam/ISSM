%Test Name: ISMIPBHO
%This test is a test from the ISMP-HOM Intercomparison project.
%Pattyn and Payne 2006
printingflag=false;

%L_list={5000.,10000.,20000.,40000.,80000.,160000.};
L_list={80000.};
results={};
minvx=[];
maxvx=[];

for i=1:length(L_list),
	L=L_list{i};
	nx=20; %numberof nodes in x direction
	ny=20;
	md=model();
	md=squaremesh(md,L,L,nx,ny);
	md=setmask(md,'',''); %ice sheet test
	md=parameterize(md,'../Par/ISMIPB.par');
	md=extrude(md,10,1.);

	md=setflowequation(md,'HO','all');

	%Create dirichlet on the bed only
	md.stressbalance.spcvx=NaN*ones(md.mesh.numberofvertices,1);
	md.stressbalance.spcvy=NaN*ones(md.mesh.numberofvertices,1);
	md.stressbalance.spcvz=NaN*ones(md.mesh.numberofvertices,1);
	pos=find(md.mesh.vertexonbase);
	md.stressbalance.spcvx(pos)=0.;
	md.stressbalance.spcvy(pos)=0.;

	%Create MPCs to have periodic boundary conditions
	posx=find(md.mesh.x==0.);
	posx2=find(md.mesh.x==max(md.mesh.x));

	posy=find(md.mesh.y==0. & md.mesh.x~=0. & md.mesh.x~=max(md.mesh.x)); %Don't take the same nodes two times
	posy2=find(md.mesh.y==max(md.mesh.y) & md.mesh.x~=0. & md.mesh.x~=max(md.mesh.x));

	md.stressbalance.vertex_pairing=[posx,posx2;posy,posy2];

	%Compute the stressbalance
	md.cluster=generic('name',oshostname(),'np',8);
	md=solve(md,'Stressbalance');

	%Plot the results and save them
	vx=(md.results.StressbalanceSolution.Vx);
	vy=(md.results.StressbalanceSolution.Vy);
	vz=(md.results.StressbalanceSolution.Vz);
	results{i}=md.results.StressbalanceSolution;
	minvx(i)=min(vx(end-md.mesh.numberofvertices2d+1:end));
	maxvx(i)=max(vx(end-md.mesh.numberofvertices2d+1:end));

	%Now plot vx, vy, vz and vx on a cross section
	plotmodel(md,'data',vx,'layer#all',md.mesh.numberoflayers,'xlim',[0 L/10^3],'ylim',[0 L/10^3],'unit','km')
	if printingflag,
		set(gcf,'Color','w')
		printmodel(['ismipbHOvx' num2str(L)],'png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
		system(['mv ismipbHOvx' num2str(L) '.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestB']);
	end
	plotmodel(md,'data',vz,'layer#all',md.mesh.numberoflayers,'xlim',[0 L/10^3],'ylim',[0 L/10^3],'unit','km')
	if printingflag,
		set(gcf,'Color','w')
		printmodel(['ismipbHOvz' num2str(L)],'png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
		system(['mv ismipbHOvz' num2str(L) '.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestB']);
	end

	if(L==5000.),
		plotmodel(md,'data',vx,'sectionvalue','../Exp/ISMIP5000.exp','layer',md.mesh.numberoflayers,...
			'resolution',[10 10],'ylim',[6 16],'xlim',[0 5000],'title','','xlabel','')
	elseif(L==10000.),
		plotmodel(md,'data',vx,'sectionvalue','../Exp/ISMIP10000.exp','layer',md.mesh.numberoflayers,...
			'resolution',[10 10],'ylim',[0 40],'xlim',[0 10000],'title','','xlabel','')
	elseif(L==20000.),
		plotmodel(md,'data',vx,'sectionvalue','../Exp/ISMIP20000.exp','layer',md.mesh.numberoflayers,...
			'resolution',[10 10],'ylim',[0 60],'xlim',[0 20000],'title','','xlabel','')
	elseif(L==40000.),
		plotmodel(md,'data',vx,'sectionvalue','../Exp/ISMIP40000.exp','layer',md.mesh.numberoflayers,...
			'resolution',[10 10],'ylim',[0 100],'xlim',[0 40000],'title','','xlabel','')
	elseif(L==80000.),
		plotmodel(md,'data',vx,'sectionvalue','../Exp/ISMIP80000.exp','layer',md.mesh.numberoflayers,...
			'resolution',[10 10],'ylim',[0 120],'xlim',[0 80000],'title','','xlabel','')
	elseif(L==160000.),
		plotmodel(md,'data',vx,'sectionvalue','../Exp/ISMIP160000.exp','layer',md.mesh.numberoflayers,...
			'resolution',[10 10],'ylim',[0 120],'xlim',[0 160000],'title','','xlabel','')
	end
	if printingflag,
		set(gcf,'Color','w')
		printmodel(['ismipbHOvxsec' num2str(L)],'png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
		system(['mv ismipbHOvxsec' num2str(L) '.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestB']);
	end
end

%Now plot the min and max values of vx for each size of the square
plot([5 10 20 40 80 160],minvx);ylim([0 14]);xlim([0 160])
if printingflag,
	set(gcf,'Color','w')
	printmodel('ismipbHOminvx','png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
	system(['mv ismipbHOminvx.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestB']);
end
plot([5 10 20 40 80 160],maxvx);ylim([0 120]);xlim([0 160])
if printingflag,
	set(gcf,'Color','w')
	printmodel('ismipbHOmaxvx','png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
	system(['mv ismipbHOmaxvx.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestB']);
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
