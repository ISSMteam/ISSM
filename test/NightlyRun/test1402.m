%Test Name: AdaptiveMeshRefinement2
%test the anisotropic mesh adaptation
printingflag=false;

%create square mesh
L=1.; %in m
nx=30; %numberof nodes in x direction
ny=30;
md=model();

%mesh adaptation loop YAMS
md=squaremesh(md,L,L,nx,ny);
u=4.*md.mesh.x-2.; v=4.*md.mesh.y-2.;
md.inversion.vel_obs=tanh(30.*(u.^2+v.^2-0.25)) ...
	+tanh(30.*((u-0.75).^2+(v-0.75).^2-0.25)) +tanh(30.*((u-0.75).^2+(v+0.75).^2-0.25)) ...
	+tanh(30.*((u+0.75).^2+(v-0.75).^2-0.25)) +tanh(30.*((u+0.75).^2+(v+0.75).^2-0.25));
plotmodel(md,'data',md.inversion.vel_obs,'data',md.inversion.vel_obs,'nlines',1,'ncols',2,'title','','figposition',[500 500 1000 500],'axis#all','equal','xlim#all',[0 1],'ylim#all',[0 1],'edgecolor#1','w'); pause(0.5);
if printingflag,
	set(gcf,'Color','w')
	printmodel('mesh2_yams1','png','margin','on','marginsize',25,'frame','off','resolution',1,'hardcopy','off');
	system(['mv mesh2_yams1.png ' ISSM_DIR '/website/doc_pdf/validation/Images/Mesh/ ']);
end

md=YamsCall(md,md.inversion.vel_obs,0.005,0.3,2.3,10.^-2);
u=4.*md.mesh.x-2.; v=4.*md.mesh.y-2.;
md.inversion.vel_obs=tanh(30.*(u.^2+v.^2-0.25)) ...
	+tanh(30.*((u-0.75).^2+(v-0.75).^2-0.25)) +tanh(30.*((u-0.75).^2+(v+0.75).^2-0.25)) ...
	+tanh(30.*((u+0.75).^2+(v-0.75).^2-0.25)) +tanh(30.*((u+0.75).^2+(v+0.75).^2-0.25));
plotmodel(md,'data',md.inversion.vel_obs,'data',md.inversion.vel_obs,'nlines',1,'ncols',2,'title','','figposition',[500 500 1000 500],'axis#all','equal','xlim#all',[0 1],'ylim#all',[0 1],'edgecolor#1','w'); pause(0.5);
if printingflag,
	set(gcf,'Color','w')
	printmodel('mesh2_yams2','png','margin','on','marginsize',25,'frame','off','resolution',1,'hardcopy','off');
	system(['mv mesh2_yams2.png ' ISSM_DIR '/website/doc_pdf/validation/Images/Mesh/ ']);
end

md=YamsCall(md,md.inversion.vel_obs,0.005,0.3,3,0.005);
u=4.*md.mesh.x-2.; v=4.*md.mesh.y-2.;
md.inversion.vel_obs=tanh(30.*(u.^2+v.^2-0.25)) ...
	+tanh(30.*((u-0.75).^2+(v-0.75).^2-0.25)) +tanh(30.*((u-0.75).^2+(v+0.75).^2-0.25)) ...
	+tanh(30.*((u+0.75).^2+(v-0.75).^2-0.25)) +tanh(30.*((u+0.75).^2+(v+0.75).^2-0.25));
plotmodel(md,'data',md.inversion.vel_obs,'data',md.inversion.vel_obs,'nlines',1,'ncols',2,'title','','figposition',[500 500 1000 500],'axis#all','equal','xlim#all',[0 1],'ylim#all',[0 1],'edgecolor#1','w'); pause(0.5);
if printingflag,
	set(gcf,'Color','w')
	printmodel('mesh2_yams3','png','margin','on','marginsize',25,'frame','off','resolution',1,'hardcopy','off');
	system(['mv mesh2_yams3.png ' ISSM_DIR '/website/doc_pdf/validation/Images/Mesh/ ']);
end
x1=md.mesh.x;
y1=md.mesh.y;

%mesh adaptation loop BAMG
md=squaremesh(md,L,L,nx,ny);
u=4.*md.mesh.x-2.; v=4.*md.mesh.y-2.;
md.inversion.vel_obs=tanh(30.*(u.^2+v.^2-0.25)) ...
	+tanh(30.*((u-0.75).^2+(v-0.75).^2-0.25)) +tanh(30.*((u-0.75).^2+(v+0.75).^2-0.25)) ...
	+tanh(30.*((u+0.75).^2+(v-0.75).^2-0.25)) +tanh(30.*((u+0.75).^2+(v+0.75).^2-0.25));
plotmodel(md,'data',md.inversion.vel_obs,'data',md.inversion.vel_obs,'nlines',1,'ncols',2,'title','','figposition',[500 500 1000 500],'axis#all','equal','xlim#all',[0 1],'ylim#all',[0 1],'edgecolor#1','w'); pause(0.5);
if printingflag,
	set(gcf,'Color','w')
	printmodel('mesh2_bamg1','png','margin','on','marginsize',25,'frame','off','resolution',1,'hardcopy','off');
	system(['mv mesh2_bamg1.png ' ISSM_DIR '/website/doc_pdf/validation/Images/Mesh/ ']);
end

md.private.bamg=NaN;
md=bamg(md,'field',md.inversion.vel_obs,'hmin',0.005,'hmax',0.3,'gradation',2.3,'err',10.^-2);
u=4.*md.mesh.x-2.; v=4.*md.mesh.y-2.;
md.inversion.vel_obs=tanh(30.*(u.^2+v.^2-0.25)) ...
	+tanh(30.*((u-0.75).^2+(v-0.75).^2-0.25)) +tanh(30.*((u-0.75).^2+(v+0.75).^2-0.25)) ...
	+tanh(30.*((u+0.75).^2+(v-0.75).^2-0.25)) +tanh(30.*((u+0.75).^2+(v+0.75).^2-0.25));
plotmodel(md,'data',md.inversion.vel_obs,'data',md.inversion.vel_obs,'nlines',1,'ncols',2,'title','','figposition',[500 500 1000 500],'axis#all','equal','xlim#all',[0 1],'ylim#all',[0 1],'edgecolor#1','w'); pause(0.5);
if printingflag,
	set(gcf,'Color','w')
	printmodel('mesh2_bamg2','png','margin','on','marginsize',25,'frame','off','resolution',1,'hardcopy','off');
	system(['mv mesh2_bamg2.png ' ISSM_DIR '/website/doc_pdf/validation/Images/Mesh/ ']);
end

md.private.bamg=NaN;
md=bamg(md,'field',md.inversion.vel_obs,'hmin',0.005,'hmax',0.3,'gradation',3,'err',0.005);
u=4.*md.mesh.x-2.; v=4.*md.mesh.y-2.;
md.inversion.vel_obs=tanh(30.*(u.^2+v.^2-0.25)) ...
	+tanh(30.*((u-0.75).^2+(v-0.75).^2-0.25)) +tanh(30.*((u-0.75).^2+(v+0.75).^2-0.25)) ...
	+tanh(30.*((u+0.75).^2+(v-0.75).^2-0.25)) +tanh(30.*((u+0.75).^2+(v+0.75).^2-0.25));
plotmodel(md,'data',md.inversion.vel_obs,'data',md.inversion.vel_obs,'nlines',1,'ncols',2,'title','','figposition',[500 500 1000 500],'axis#all','equal','xlim#all',[0 1],'ylim#all',[0 1],'edgecolor#1','w'); pause(0.5);
if printingflag,
	set(gcf,'Color','w')
	printmodel('mesh2_bamg3','png','margin','on','marginsize',25,'frame','off','resolution',1,'hardcopy','off');
	system(['mv mesh2_bamg3.png ' ISSM_DIR '/website/doc_pdf/validation/Images/Mesh/ ']);
end

md.private.bamg=NaN;
md=bamg(md,'field',md.inversion.vel_obs,'hmin',0.005,'hmax',0.3,'gradation',1.5,'err',0.003,'anisomax',1);
u=4.*md.mesh.x-2.; v=4.*md.mesh.y-2.;
md.inversion.vel_obs=tanh(30.*(u.^2+v.^2-0.25)) ...
	+tanh(30.*((u-0.75).^2+(v-0.75).^2-0.25)) +tanh(30.*((u-0.75).^2+(v+0.75).^2-0.25)) ...
	+tanh(30.*((u+0.75).^2+(v-0.75).^2-0.25)) +tanh(30.*((u+0.75).^2+(v+0.75).^2-0.25));
plotmodel(md,'data',md.inversion.vel_obs,'data',md.inversion.vel_obs,'nlines',1,'ncols',2,'title','','figposition',[500 500 1000 500],'axis#all','equal','xlim#all',[0 1],'ylim#all',[0 1],'edgecolor#1','w'); pause(0.5);
if printingflag,
	set(gcf,'Color','w')
	printmodel('mesh2_bamgiso','png','margin','on','marginsize',25,'frame','off','resolution',1,'hardcopy','off');
	system(['mv mesh2_bamgiso.png ' ISSM_DIR '/website/doc_pdf/validation/Images/Mesh/ ']);
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
