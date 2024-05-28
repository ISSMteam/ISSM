%Test Name: ISMIPE
%This test is a test from the ISMP-HOM Intercomparison project.
%TestE 
%Four tests to run: - Pattyn frozen
%                   - Stokes frozen
%                   - Pattyn with some sliding
%                   - Stokes with some sliding
printingflag=false;
results={};

for i=1:4,
	Lx=10.; %in m
	Ly=5000.; %in m
	nx=3; %number of nodes in x direction
	ny=51;
	md=model();
	md=squaremesh(md,Lx,Ly,nx,ny);
	md=setmask(md,'',''); %ice sheet test
	md=parameterize(md,'../Par/ISMIPE.par');
	md=extrude(md,10,1.);

	if i==1 | i==3,
		md=setflowequation(md,'HO','all');
	elseif i==2 | i==4,
		md=setflowequation(md,'FS','all');
	end

	%Create MPCs to have periodic boundary conditions
	posx=find(md.mesh.x==0.);
	posx2=find(md.mesh.x==max(md.mesh.x));
	md.stressbalance.vertex_pairing=[posx,posx2];

	%Create spcs on the bed 
	pos=find(md.mesh.vertexonbase);
	md.stressbalance.spcvx=NaN*ones(md.mesh.numberofvertices,1);
	md.stressbalance.spcvy=NaN*ones(md.mesh.numberofvertices,1);
	md.stressbalance.spcvz=NaN*ones(md.mesh.numberofvertices,1);
	md.stressbalance.spcvx(pos)=0.;
	md.stressbalance.spcvy(pos)=0.;
	md.stressbalance.spcvz(pos)=0.;

	%Remove the spc where there is some sliding (case 3 and 4):
	if i==3 | i==4,
		pos=find(md.mesh.y/max(md.mesh.y)>=0.44 & md.mesh.y/max(md.mesh.y)<=0.5);
		md.stressbalance.spcvx(pos)=NaN;
		md.stressbalance.spcvy(pos)=NaN;
		md.stressbalance.spcvz(pos)=NaN;
	end

	%Compute the stressbalance
	md.cluster=generic('name',oshostname(),'np',8);
	md=solve(md,'Stressbalance');

	vx=(md.results.StressbalanceSolution.Vx);
	vy=(md.results.StressbalanceSolution.Vy);
	vz=(md.results.StressbalanceSolution.Vz);
	results{i}=md.results.StressbalanceSolution;

	if i==1,
		plotmodel(md,'data',vy,'ylim',[-10 80],'layer',md.mesh.numberoflayers,'sectionvalue','../Exp/ISMIPE.exp','resolution',[10 10],'title','','xlabel','')
		if printingflag,
			set(gcf,'Color','w')
			printmodel('ismipeHOvxfrozen','png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
			system(['mv ismipeHOvxfrozen.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestE']);
		end
	elseif i==2,
		plotmodel(md,'data',vy,'ylim',[-10 80],'layer',md.mesh.numberoflayers,'sectionvalue','../Exp/ISMIPE.exp','resolution',[10 10],'title','','xlabel','')
		if printingflag,
			set(gcf,'Color','w')
			printmodel('ismipeFSvxfrozen','png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
			system(['mv ismipeFSvxfrozen.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestE']);
		end
	elseif i==3,
		plotmodel(md,'data',vy,'ylim',[-50 200],'layer',md.mesh.numberoflayers,'sectionvalue','../Exp/ISMIPE.exp','resolution',[10 10],'title','','xlabel','')
		if printingflag,
			set(gcf,'Color','w')
			printmodel('ismipeHOvxsliding','png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
			system(['mv ismipeHOvxsliding.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestE']);
		end
	elseif i==4,
		plotmodel(md,'data',vy,'ylim',[-50 200],'layer',md.mesh.numberoflayers,'sectionvalue','../Exp/ISMIPE.exp','resolution',[10 10],'title','','xlabel','')
		if printingflag,
			set(gcf,'Color','w')
			printmodel('ismipeFSvxsliding','png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
			system(['mv ismipeFSvxsliding.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestE']);
		end
	end
end

%Fields and tolerances to track changes
field_names     ={...
	'VyPattynSliding','VzPattynSliding',...
	'VxStokesSliding','VyStokesSliding','VzStokesSliding',...
	'VyPattynFrozen','VzPattynFrozen',...
	'VxStokesFrozen','VyStokesFrozen','VzStokesFrozen'
};
field_tolerances={...
	1e-05,1e-05,...
	1e-05,1e-06,1e-06,...
	1e-05,1e-04,...
	1e-05,1e-05,1e-06,...
};
field_values={};
for i=1:4,
	result=results{i};
	field_values={field_values{:},...
		(result.Vx),...
		(result.Vy),...
		(result.Vz),...
		};
end
