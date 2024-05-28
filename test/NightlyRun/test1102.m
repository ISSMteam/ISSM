%Test Name: ISMIPAFS
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

%	%Find elements at the corner and extract model
%	posnodes=find((md.mesh.x==0. | md.mesh.x==max(md.mesh.x)) & (md.mesh.y==0. | md.mesh.y==max(md.mesh.y)));
%	[a,b]=find(ismember(md.mesh.elements,posnodes));
%	elements=ones(md.mesh.numberofelements,1);
%	elements(a)=0;
%	md=modelextract(md,elements);

	md=parameterize(md,'../Par/ISMIPA.par');
	md=extrude(md,10,1.);
	md=setflowequation(md,'FS','all');

	%Create dirichlet on the bed only
	pos=find(md.mesh.vertexonbase);
	md.stressbalance.spcvx(pos)=0.;
	md.stressbalance.spcvy(pos)=0.;
	md.stressbalance.spcvz(pos)=0.;

%	%Create MPCs to have periodic boundary conditions
%	posx=find(md.mesh.x==0.);
%	posx2=find(md.mesh.x==max(md.mesh.x));
%	posx=find(md.mesh.x==0. & md.mesh.y~=0. & md.mesh.y~=max(md.mesh.y) & ~md.mesh.vertexonbase);
%	posx2=find(md.mesh.x==max(md.mesh.x) &  md.mesh.y~=0. & md.mesh.y~=max(md.mesh.y) & ~md.mesh.vertexonbase);

%	posy=find(md.mesh.y==0. & md.mesh.x~=0. & md.mesh.x~=max(md.mesh.x) & ~md.mesh.vertexonbase); %Don't take the same nodes two times
%	posy2=find(md.mesh.y==max(md.mesh.y) & md.mesh.x~=0. & md.mesh.x~=max(md.mesh.x) & ~md.mesh.vertexonbase);

%	md.stressbalance.vertex_pairing=[posx,posx2;posy,posy2];

	%Compute the stressbalance
	md.stressbalance.abstol=NaN;
	md.stressbalance.reltol=NaN;
	md.stressbalance.restol=1.;
	md.cluster=generic('name',oshostname(),'np',8);
	md=solve(md,'Stressbalance');

	%Plot the results and save them
	vx=(md.results.StressbalanceSolution.Vx);
	vy=(md.results.StressbalanceSolution.Vy);
	vz=(md.results.StressbalanceSolution.Vz);
	pressure=(md.results.StressbalanceSolution.Pressure);
	results{i}=md.results.StressbalanceSolution;
	minvx(i)=min(vx(end-md.mesh.numberofvertices2d+1:end));
	maxvx(i)=max(vx(end-md.mesh.numberofvertices2d+1:end));

	%Now plot vx, vy, vz and vx on a cross section
	plotmodel(md,'data',vx,'layer#all',md.mesh.numberoflayers,'xlim',[0 L/10^3],'ylim',[0 L/10^3],'unit','km','figure',2)
	if printingflag,
		set(gcf,'Color','w')
		printmodel(['ismipaFSvx' num2str(L)],'png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
		system(['mv ismipaFSvx' num2str(L) '.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestA']);
	end
	plotmodel(md,'data',vy,'layer#all',md.mesh.numberoflayers,'xlim',[0 L/10^3],'ylim',[0 L/10^3],'unit','km','figure',3)
	if printingflag,
		set(gcf,'Color','w')
		printmodel(['ismipaFSvy' num2str(L)],'png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
		system(['mv ismipaFSvy' num2str(L) '.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestA']);
	end
	plotmodel(md,'data',vz,'layer#all',md.mesh.numberoflayers,'xlim',[0 L/10^3],'ylim',[0 L/10^3],'unit','km','figure',4)
	if printingflag,
		set(gcf,'Color','w')
		printmodel(['ismipaFSvz' num2str(L)],'png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
		system(['mv ismipaFSvz' num2str(L) '.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestA']);
	end

	if(L==5000.),
		plotmodel(md,'data',vx,'sectionvalue','../Exp/ISMIP5000.exp','layer',md.mesh.numberoflayers,...
			'resolution',[10 10],'ylim',[10 18],'xlim',[0 5000],'title','','xlabel','')
	elseif(L==10000.),
		plotmodel(md,'data',vx,'sectionvalue','../Exp/ISMIP10000.exp','layer',md.mesh.numberoflayers,...
			'resolution',[10 10],'ylim',[10 30],'xlim',[0 10000],'title','','xlabel','')
	elseif(L==20000.),
		plotmodel(md,'data',vx,'sectionvalue','../Exp/ISMIP20000.exp','layer',md.mesh.numberoflayers,...
			'resolution',[10 10],'ylim',[0 50],'xlim',[0 20000],'title','','xlabel','')
	elseif(L==40000.),
		plotmodel(md,'data',vx,'sectionvalue','../Exp/ISMIP40000.exp','layer',md.mesh.numberoflayers,...
			'resolution',[10 10],'ylim',[0 80],'xlim',[0 40000],'title','','xlabel','')
	elseif(L==80000.),
		plotmodel(md,'data',vx,'sectionvalue','../Exp/ISMIP80000.exp','layer',md.mesh.numberoflayers,...
			'resolution',[10 10],'ylim',[0 100],'xlim',[0 80000],'title','','xlabel','')
	elseif(L==160000.),
		plotmodel(md,'data',vx,'sectionvalue','../Exp/ISMIP160000.exp','layer',md.mesh.numberoflayers,...
			'resolution',[10 10],'ylim',[0 120],'xlim',[0 160000],'title','','xlabel','')
	end
	if printingflag,
		set(gcf,'Color','w')
		printmodel(['ismipaFSvxsec' num2str(L)],'png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
		system(['mv ismipaFSvxsec' num2str(L) '.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestA']);
	end
end

%Now plot the min and max values of vx for each size of the square
plot([5 10 20 40 80 160],minvx);ylim([0 18])
if printingflag,
	set(gcf,'Color','w')
	printmodel('ismipaFSminvx','png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
	system(['mv ismipaFSminvx.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestA']);
end
plot([5 10 20 40 80 160],maxvx);ylim([0 120])
if printingflag,
	set(gcf,'Color','w')
	printmodel('ismipaFSmaxvx','png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
	system(['mv ismipaFSmaxvx.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestA']);
end

%Fields and tolerances to track changes
field_names     ={...
	'Vx80km','Vy80km','Vz80km'
};
field_tolerances={...
	1e-12,1e-11,1e-12,...
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
