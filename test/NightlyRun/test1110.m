%Test Name: ISMIPF
%This test is a test from the ISMP-HOM Intercomparison project.
%TestF
printingflag=false;
results={};

for i=3,  %1:4,
	L=100000.; %in m
	nx=30; %numberof nodes in x direction
	ny=30;
	md=model();
	md=squaremesh(md,L,L,nx,ny);
%	md=triangle(md,'../Exp/SquareISMIP.exp',5500.);
	md=setmask(md,'',''); %ice sheet test
	md=parameterize(md,'../Par/ISMIPF.par');
	md=extrude(md,4,1.);

	if (i==1 | i==2),
		md=setflowequation(md,'HO','all');
	else
		md=setflowequation(md,'FS','all');
	end

	md.stressbalance.spcvx=NaN*ones(md.mesh.numberofvertices,1);
	md.stressbalance.spcvy=NaN*ones(md.mesh.numberofvertices,1);
	md.stressbalance.spcvz=NaN*ones(md.mesh.numberofvertices,1);
	if (i==1 | i==3),
		%Create dirichlet on the bed if no slip
		pos=find(md.mesh.vertexonbase);
		md.stressbalance.spcvx(pos)=0.;
		md.stressbalance.spcvy(pos)=0.;
		md.stressbalance.spcvz(pos)=0.;
	else
		pos=find(md.mesh.vertexonbase & (md.mesh.x==0. | md.mesh.x==max(md.mesh.x)) & (md.mesh.y==0. | md.mesh.y==max(md.mesh.y)));
		md.stressbalance.spcvx(pos)=100.; %because we need a dirichlet somewhere
		md.stressbalance.spcvy(pos)=0.;
		md.stressbalance.spcvz(pos)=0.;
	end
	pos=find(~md.mesh.vertexonbase);
	md.thermal.spctemperature(pos)=255.;

	%Create MPCs to have periodic boundary conditions
	posx=find(md.mesh.x==0.);
	posx2=find(md.mesh.x==max(md.mesh.x));

	posy=find(md.mesh.y==0.);
	posy2=find(md.mesh.y==max(md.mesh.y));

	md.stressbalance.vertex_pairing=[posx,posx2;posy,posy2];
	md.masstransport.vertex_pairing=[posx,posx2;posy,posy2];

	md.timestepping.time_step=3.;
	md.timestepping.final_time=300.;
	md.settings.output_frequency=50;
	md.masstransport.stabilization=1;
	md.stressbalance.maxiter=1;

	%Compute the stressbalance
	md.cluster=generic('name',oshostname(),'np',8);
	md.verbose=verbose('convergence',true,'solution',true);
	md=solve(md,'Transient');

	%save the results
	results{i}=md.results.TransientSolution(end);

	%Now plot vx and delta surface
	if (i==1 | i==3),
		plotmodel(md,'data',(md.results.TransientSolution(end).Vx),'layer',md.mesh.numberoflayers,'sectionvalue','../Exp/ISMIP100000.exp','title','','xlabel','','ylabel','Velocity (m/yr)','linewidth',3,'grid','on','unit','km','ylim',[91 100])
	elseif (i==2 | i==4),
		plotmodel(md,'data',(md.results.TransientSolution(end).Vx),'layer',md.mesh.numberoflayers,'sectionvalue','../Exp/ISMIP100000.exp','title','','xlabel','','ylabel','Velocity (m/yr)','linewidth',3,'grid','on','unit','km','ylim',[185 200])
	end
	if printingflag,
		set(gcf,'Color','w')
		if i==1,
			printmodel('ismipfHOvxfrozen','png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
			system(['mv ismipfHOvxfrozen.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestF']);
		elseif i==2,
			printmodel('ismipfHOvxsliding','png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
			system(['mv ismipfHOvxsliding.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestF']);
		elseif i==3,
			printmodel('ismipfFSvxfrozen','png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
			system(['mv ismipfFSvxfrozen.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestF']);
		elseif i==4,
			printmodel('ismipfFSvxsliding','png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
			system(['mv ismipfFSvxsliding.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestF']);
		end
	end

	plotmodel(md,'data',(md.results.TransientSolution(end).Surface)-md.geometry.surface,'layer',md.mesh.numberoflayers,'sectionvalue','../Exp/ISMIP100000.exp','title','','xlabel','','ylabel','Surface (m)','linewidth',3,'grid','on','unit','km','ylim',[-30 50])
	if printingflag,
		set(gcf,'Color','w')
		if i==1,
			printmodel('ismipfHOdeltasurfacefrozen','png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
			system(['mv ismipfHOdeltasurfacefrozen.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestF']);
		elseif i==2,
			printmodel('ismipfHOdeltasurfacesliding','png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
			system(['mv ismipfHOdeltasurfacesliding.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestF']);
		elseif i==3,
			printmodel('ismipfFSdeltasurfacefrozen','png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
			system(['mv ismipfFSdeltasurfacefrozen.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestF']);
		elseif i==4,
			printmodel('ismipfFSdeltasurfacesliding','png','margin','on','marginsize',25,'frame','off','resolution',1.5,'hardcopy','off');
			system(['mv ismipfFSdeltasurfacesliding.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestF']);
		end
	end
end

%Fields and tolerances to track changes
field_names     ={...
	'VxPattynFrozen','VyPattynFrozen','VzPattynFrozen','SurfacePattynFrozen',...
	'VxPattynSliding','VyPattynSliding','VzPattynSliding','SurfacePattynSliding',...
	'VxStokesFrozen','VyStokesFrozen','VzStokesFrozen','SurfaceStokesFrozen',...
	'VxStokesSliding','VyStokesSliding','VzStokesSliding','SurfaceStokesSliding'
};
field_tolerances={...
	1e-10,1e-09,1e-09,1e-10,...
	1e-10,1e-09,1e-09,1e-10,...
	1e-08,1e-09,1e-08,1e-09,...
	1e-08,2e-09,1e-08,1e-09,...
};
field_values={};
for i=1:4,
	result=results{i};
	field_values={field_values{:},...
		(result.Vx),...
		(result.Vy),...
		(result.Vz),...
		(result.Surface)-md.geometry.surface,...
		};
end
