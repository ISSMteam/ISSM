%Test Name: ISMIPDHO
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
	nx=30; %numberof nodes in x direction
	ny=30;
	md=model();
	md=squaremesh(md,L,L,nx,ny);
	md=setmask(md,'',''); %ice sheet test
	md=parameterize(md,'../Par/ISMIPD.par');
	md=extrude(md,10,1.);

	md=setflowequation(md,'HO','all');

	%We need one grd on dirichlet: the 4 corners are set to zero
	md.stressbalance.spcvx=NaN*ones(md.mesh.numberofvertices,1);
	md.stressbalance.spcvy=NaN*ones(md.mesh.numberofvertices,1);
	md.stressbalance.spcvz=NaN*ones(md.mesh.numberofvertices,1);

	%Create MPCs to have periodic boundary conditions
	posx=find(md.mesh.x==0. & ~(md.mesh.y==0. & md.mesh.vertexonbase) & ~(md.mesh.y==L & md.mesh.vertexonbase));
	posx2=find(md.mesh.x==max(md.mesh.x) & ~(md.mesh.y==0. & md.mesh.vertexonbase) & ~(md.mesh.y==L & md.mesh.vertexonbase));

	posy=find(md.mesh.y==0. & md.mesh.x~=0. & md.mesh.x~=max(md.mesh.x)); %Don't take the same nodes two times
	posy2=find(md.mesh.y==max(md.mesh.y) & md.mesh.x~=0. & md.mesh.x~=max(md.mesh.x));

	md.stressbalance.vertex_pairing=[posx,posx2;posy,posy2];

	%Add spc on the corners
	pos=find((md.mesh.x==0. | md.mesh.x==L) & (md.mesh.y==0. | md.mesh.y==L) & md.mesh.vertexonbase);
	md.stressbalance.spcvy(:)=0.;
	md.stressbalance.spcvx(pos)=0.;
	if(L==5000.),
		md.stressbalance.spcvx(pos)=16.0912;
	elseif(L==10000.),
		md.stressbalance.spcvx(pos)=16.52;
	elseif(L==20000.),
		md.stressbalance.spcvx(pos)=17.77;
	elseif(L==40000.),
		md.stressbalance.spcvx(pos)=19.88;
	elseif(L==80000.),
		md.stressbalance.spcvx(pos)=18.65;
	elseif(L==160000.),
		md.stressbalance.spcvx(pos)=16.91;
	end
	
	%Spc the bed at zero for vz
	pos=find(md.mesh.vertexonbase);
	md.stressbalance.spcvz(pos)=0.;

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
	plotmodel(md,'data',vx,'layer#all',md.mesh.numberoflayers,'xlim',[0 L/10^3],'ylim',[0 L/10^3],'unit','km','figure',2)
	if printingflag,
		set(gcf,'Color','w')
		printmodel(['ismipdHOvx' num2str(L)],'png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
		system(['mv ismipdHOvx' num2str(L) '.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestD']);
	end
	plotmodel(md,'data',vz,'layer#all',md.mesh.numberoflayers,'xlim',[0 L/10^3],'ylim',[0 L/10^3],'unit','km','figure',3)
	if printingflag,
		set(gcf,'Color','w')
		printmodel(['ismipdHOvz' num2str(L)],'png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
		system(['mv ismipdHOvz' num2str(L) '.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestD']);
	end

	if(L==5000.),
		plotmodel(md,'data',vx,'sectionvalue','../Exp/ISMIP5000.exp','layer',md.mesh.numberoflayers,...
			'resolution',[10 10],'ylim',[0 20],'xlim',[0 5000],'title','','xlabel','','figure',4)
	elseif(L==10000.),
		plotmodel(md,'data',vx,'sectionvalue','../Exp/ISMIP10000.exp','layer',md.mesh.numberoflayers,...
			'resolution',[10 10],'ylim',[0 20],'xlim',[0 10000],'title','','xlabel','','figure',4)
	elseif(L==20000.),
		plotmodel(md,'data',vx,'sectionvalue','../Exp/ISMIP20000.exp','layer',md.mesh.numberoflayers,...
			'resolution',[10 10],'ylim',[0 30],'xlim',[0 20000],'title','','xlabel','','figure',4)
	elseif(L==40000.),
		plotmodel(md,'data',vx,'sectionvalue','../Exp/ISMIP40000.exp','layer',md.mesh.numberoflayers,...
			'resolution',[10 10],'ylim',[10 60],'xlim',[0 40000],'title','','xlabel','','figure',4)
	elseif(L==80000.),
		plotmodel(md,'data',vx,'sectionvalue','../Exp/ISMIP80000.exp','layer',md.mesh.numberoflayers,...
			'resolution',[10 10],'ylim',[0 200],'xlim',[0 80000],'title','','xlabel','','figure',4)
	elseif(L==160000.),
		plotmodel(md,'data',vx,'sectionvalue','../Exp/ISMIP160000.exp','layer',md.mesh.numberoflayers,...
			'resolution',[10 10],'ylim',[0 400],'xlim',[0 160000],'title','','xlabel','','figure',4)
	end
	if printingflag,
		set(gcf,'Color','w')
		printmodel(['ismipdHOvxsec' num2str(L)],'png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
		system(['mv ismipdHOvxsec' num2str(L) '.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestD']);
	end
end

%Now plot the min and max values of vx for each size of the square
plot([5 10 20 40 80 160],minvx);ylim([2 18]);xlim([0 160])
if printingflag,
	set(gcf,'Color','w')
	printmodel('ismipdHOminvx','png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
	system(['mv ismipdHOminvx.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestD']);
end
plot([5 10 20 40 80 160],maxvx);ylim([0 300]);xlim([0 160])
if printingflag,
	set(gcf,'Color','w')
	printmodel('ismipdHOmaxvx','png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
	system(['mv ismipdHOmaxvx.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestD']);
end

%Fields and tolerances to track changes
field_names     ={...
	'Vx80km','Vy80km','Vz80km'
};
field_tolerances={...
	1e-08,1e-08,1e-07,...
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
