%Test Name: AdaptiveMeshRefinement1
%test the anisotropic mesh adaptation
%function to capture = exp(-(sqrt((md.mesh.x+0.1).^2+(md.mesh.y+0.1).^2)-0.75).^2*10.^6)+((md.mesh.x+0.1).^2+(md.mesh.y+0.1).^2)/2.;
printingflag=false;

%create square mesh
L=1.; %in m
nx=70; %numberof nodes in x direction
ny=70;
md=model();

%mesh adaptation loop YAMS
md=squaremesh(md,L,L,nx,ny);
md.inversion.vel_obs=exp(-(sqrt((md.mesh.x+0.1).^2+(md.mesh.y+0.1).^2)-0.75).^2*10.^6)+((md.mesh.x+0.1).^2+(md.mesh.y+0.1).^2)/2.;
plotmodel(md,'data',md.inversion.vel_obs,'data',md.inversion.vel_obs,'nlines',1,'ncols',2,'title','','figposition',[500 500 1000 500],'axis#all','equal','xlim#all',[0 1],'ylim#all',[0 1],'edgecolor#1','w'); pause(0.5);
if printingflag,
	set(gcf,'Color','w')
	printmodel('mesh1_yams1','png','margin','on','marginsize',25,'frame','off','resolution',1,'hardcopy','off');
	system(['mv mesh1_yams1.png ' ISSM_DIR '/website/doc_pdf/validation/Images/Mesh/ ']);
end

md=YamsCall(md,md.inversion.vel_obs,0.001,0.3,1.3,10.^-4);
md.inversion.vel_obs=exp(-(sqrt((md.mesh.x+0.1).^2+(md.mesh.y+0.1).^2)-0.75).^2*10.^6)+((md.mesh.x+0.1).^2+(md.mesh.y+0.1).^2)/2.;
plotmodel(md,'data',md.inversion.vel_obs,'data',md.inversion.vel_obs,'nlines',1,'ncols',2,'title','','figposition',[500 500 1000 500],'axis#all','equal','xlim#all',[0 1],'ylim#all',[0 1],'edgecolor#1','w'); pause(0.5);
if printingflag,
	set(gcf,'Color','w')
	printmodel('mesh1_yams2','png','margin','on','marginsize',25,'frame','off','resolution',1,'hardcopy','off');
	system(['mv mesh1_yams2.png ' ISSM_DIR '/website/doc_pdf/validation/Images/Mesh/ ']);
end

md=YamsCall(md,md.inversion.vel_obs,0.001,0.3,2.5,0.008);
md.inversion.vel_obs=exp(-(sqrt((md.mesh.x+0.1).^2+(md.mesh.y+0.1).^2)-0.75).^2*10.^6)+((md.mesh.x+0.1).^2+(md.mesh.y+0.1).^2)/2.;
plotmodel(md,'data',md.inversion.vel_obs,'data',md.inversion.vel_obs,'nlines',1,'ncols',2,'title','','figposition',[500 500 1000 500],'axis#all','equal','xlim#all',[0 1],'ylim#all',[0 1],'edgecolor#1','w'); pause(0.5);
if printingflag,
	set(gcf,'Color','w')
	printmodel('mesh1_yams3','png','margin','on','marginsize',25,'frame','off','resolution',1,'hardcopy','off');
	system(['mv mesh1_yams3.png ' ISSM_DIR '/website/doc_pdf/validation/Images/Mesh/ ']);
end
x1=md.mesh.x;
y1=md.mesh.y;

%mesh adaptation loop BAMG
md=squaremesh(md,L,L,nx,ny);
md.inversion.vel_obs=exp(-(sqrt((md.mesh.x+0.1).^2+(md.mesh.y+0.1).^2)-0.75).^2*10.^6)+((md.mesh.x+0.1).^2+(md.mesh.y+0.1).^2)/2.;
plotmodel(md,'data',md.inversion.vel_obs,'data',md.inversion.vel_obs,'nlines',1,'ncols',2,'title','','figposition',[500 500 1000 500],'axis#all','equal','xlim#all',[0 1],'ylim#all',[0 1],'edgecolor#1','w'); pause(0.5);
if printingflag,
	set(gcf,'Color','w')
	printmodel('mesh1_bamg1','png','margin','on','marginsize',25,'frame','off','resolution',1,'hardcopy','off');
	system(['mv mesh1_bamg1.png ' ISSM_DIR '/website/doc_pdf/validation/Images/Mesh/ ']);
end

md.private.bamg=NaN;
md=bamg(md,'field',md.inversion.vel_obs,'hmin',0.001,'hmax',0.3,'gradation',1.3,'err',10.^-4);
md.inversion.vel_obs=exp(-(sqrt((md.mesh.x+0.1).^2+(md.mesh.y+0.1).^2)-0.75).^2*10.^6)+((md.mesh.x+0.1).^2+(md.mesh.y+0.1).^2)/2.;
plotmodel(md,'data',md.inversion.vel_obs,'data',md.inversion.vel_obs,'nlines',1,'ncols',2,'title','','figposition',[500 500 1000 500],'axis#all','equal','xlim#all',[0 1],'ylim#all',[0 1],'edgecolor#1','w'); pause(0.5);
if printingflag,
	set(gcf,'Color','w')
	printmodel('mesh1_bamg2','png','margin','on','marginsize',25,'frame','off','resolution',1,'hardcopy','off');
	system(['mv mesh1_bamg2.png ' ISSM_DIR '/website/doc_pdf/validation/Images/Mesh/ ']);
end

md.private.bamg=NaN;
md=bamg(md,'field',md.inversion.vel_obs,'hmin',0.001,'hmax',0.3,'gradation',2.5,'err',0.008);
md.inversion.vel_obs=exp(-(sqrt((md.mesh.x+0.1).^2+(md.mesh.y+0.1).^2)-0.75).^2*10.^6)+((md.mesh.x+0.1).^2+(md.mesh.y+0.1).^2)/2.;
plotmodel(md,'data',md.inversion.vel_obs,'data',md.inversion.vel_obs,'nlines',1,'ncols',2,'title','','figposition',[500 500 1000 500],'axis#all','equal','xlim#all',[0 1],'ylim#all',[0 1],'edgecolor#1','w'); pause(0.5);
if printingflag,
	set(gcf,'Color','w')
	printmodel('mesh1_bamg3','png','margin','on','marginsize',25,'frame','off','resolution',1,'hardcopy','off');
	system(['mv mesh1_bamg3.png ' ISSM_DIR '/website/doc_pdf/validation/Images/Mesh/ ']);
end
x2=md.mesh.x;
y2=md.mesh.y;

%Fields and tolerances to track changes
field_names     ={'xyams','yyams','xbamg','ybamg'};
field_tolerances={1e-13,1e-13,1e-13,1e-13};
field_values={...
	x1, y1,...
	x2, y2,...
	};
