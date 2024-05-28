%Test Name: ISMIPCHO
%This test is a test from the ISMP-HOM Intercomparison project.
%Pattyn and Payne 2006
printingflag=false;

%L_list={5000.,10000.,20000.,40000.,80000.,160000.};
L_list={80000.};
results={};
minvx=[];
maxvx=[];

for i=1:length(L_list),
	L=L_list{i};  %in m (3 times the desired length for BC problems)  
	nx=30; %number of nodes in x direction
	ny=30;
	md=model();
	md=squaremesh(md,L,L,nx,ny);
	md=setmask(md,'',''); %ice sheet test
	md=parameterize(md,'../Par/ISMIPC.par');
	md=extrude(md,10,1.);

	md=setflowequation(md,'HO','all');

	%Create MPCs to have periodic boundary conditions
	md.stressbalance.spcvx=NaN*ones(md.mesh.numberofvertices,1);
	md.stressbalance.spcvy=NaN*ones(md.mesh.numberofvertices,1);
	md.stressbalance.spcvz=NaN*ones(md.mesh.numberofvertices,1);

	posx=find(md.mesh.x==0. & md.mesh.y~=0. & md.mesh.y~=L);
	posx2=find(md.mesh.x==L & md.mesh.y~=0. & md.mesh.y~=L);

	posy=find(md.mesh.y==0. & md.mesh.x~=0. & md.mesh.x~=L); %Don't take the same nodes two times
	posy2=find(md.mesh.y==L & md.mesh.x~=0. & md.mesh.x~=L);

	md.stressbalance.vertex_pairing=[posx,posx2;posy,posy2];

	%Add spc on the corners
	pos=find((md.mesh.x==0. | md.mesh.x==L) & (md.mesh.y==0. | md.mesh.y==L) & md.mesh.vertexonbase);
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
		printmodel(['ismipcHOvx' num2str(L)],'png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
		system(['mv ismipcHOvx' num2str(L) '.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestC']);
	end
	plotmodel(md,'data',vy,'layer#all',md.mesh.numberoflayers,'xlim',[0 L/10^3],'ylim',[0 L/10^3],'unit','km','figure',3)
	if printingflag,
		set(gcf,'Color','w')
		printmodel(['ismipcHOvy' num2str(L)],'png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
		system(['mv ismipcHOvy' num2str(L) '.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestC']);
	end
	plotmodel(md,'data',vz,'layer#all',md.mesh.numberoflayers,'xlim',[0 L/10^3],'ylim',[0 L/10^3],'unit','km','figure',4)
	if printingflag,
		set(gcf,'Color','w')
		printmodel(['ismipcHOvz' num2str(L)],'png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
		system(['mv ismipcHOvz' num2str(L) '.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestC']);
	end

	if(L==5000.),
		plotmodel(md,'data',vx,'sectionvalue','../Exp/ISMIP5000.exp','layer',md.mesh.numberoflayers,...
			'resolution',[10 10],'ylim',[0 20],'xlim',[0 5000],'title','','xlabel','','figure',5)
	elseif(L==10000.),
		plotmodel(md,'data',vx,'sectionvalue','../Exp/ISMIP10000.exp','layer',md.mesh.numberoflayers,...
			'resolution',[10 10],'ylim',[13 18],'xlim',[0 10000],'title','','xlabel','')
	elseif(L==20000.),
		plotmodel(md,'data',vx,'sectionvalue','../Exp/ISMIP20000.exp','layer',md.mesh.numberoflayers,...
			'resolution',[10 10],'ylim',[14 22],'xlim',[0 20000],'title','','xlabel','')
	elseif(L==40000.),
		plotmodel(md,'data',vx,'sectionvalue','../Exp/ISMIP40000.exp','layer',md.mesh.numberoflayers,...
			'resolution',[10 10],'ylim',[10 40],'xlim',[0 40000],'title','','xlabel','')
	elseif(L==80000.),
		plotmodel(md,'data',vx,'sectionvalue','../Exp/ISMIP80000.exp','layer',md.mesh.numberoflayers,...
			'resolution',[10 10],'ylim',[0 80],'xlim',[0 80000],'title','','xlabel','')
	elseif(L==160000.),
		plotmodel(md,'data',vx,'sectionvalue','../Exp/ISMIP160000.exp','layer',md.mesh.numberoflayers,...
			'resolution',[10 10],'ylim',[0 200],'xlim',[0 160000],'title','','xlabel','')
	end
	if printingflag,
		set(gcf,'Color','w')
		printmodel(['ismipcHOvxsec' num2str(L)],'png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
		system(['mv ismipcHOvxsec' num2str(L) '.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestC']);
	end
end

%Now plot the min and max values of vx for each size of the square
plot([5 10 20 40 80 160],minvx);ylim([4 18]);xlim([0 160])
if printingflag,
	set(gcf,'Color','w')
	printmodel('ismipcHOminvx','png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
	system(['mv ismipcHOminvx.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestC']);
end
plot([5 10 20 40 80 160],maxvx);ylim([0 200]); xlim([0 160])
if printingflag,
	set(gcf,'Color','w')
	printmodel('ismipcHOmaxvx','png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
	system(['mv ismipcHOmaxvx.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestC']);
end

%Fields and tolerances to track changes
field_names     ={...
	'Vx80km','Vy80km','Vz80km'
};
field_tolerances={...
	1e-09,1e-08,1e-08,...
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
