function plot_referential(md,options,width,i,data)
%PLOT_PRESSURELOAD - plot segment on neumann BC
%
%   Usage:
%      plot_referential(md,options,width,i);
%
%   See also: PLOTMODEL

%plot mesh boundaries
subplot(width,width,i); 

%process mesh and data
[x y z elements is2d isplanet]=processmesh(md,[],options);
referential=md.stressbalance.referential;

Xhat=md.stressbalance.referential(:,1:3);
pos=find(sum(isnan(Xhat),2));
Xhat(pos,:)=repmat([1 0 0],size(pos,1),1);
Xhatnorm=sqrt(Xhat(:,1).^2+Xhat(:,2).^2+Xhat(:,3).^2);
Xhat=Xhat./[Xhatnorm Xhatnorm Xhatnorm];

Zhat=md.stressbalance.referential(:,4:6);
pos=find(sum(isnan(Zhat),2));
Zhat(pos,:)=repmat([0 0 1],size(pos,1),1);
Zhatnorm=sqrt(Zhat(:,1).^2+Zhat(:,2).^2+Zhat(:,3).^2);
Zhat=Zhat./[Zhatnorm Zhatnorm Zhatnorm];

Yhat=cross(Zhat,Xhat);

if dimension(md.mesh)==2,

	%plot mesh
	A=elements(:,1); B=elements(:,2); C=elements(:,3); 
	h1=patch( 'Faces', [A B C], 'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
	hold on;

	xstart=x;
	ystart=y;
	zstart=z;
	edgex=max(md.mesh.x(elements),[],2)-min(md.mesh.x(elements),[],2);
	len=min(edgex)/1.5;
	%plot X axis
	xend=xstart+len*Xhat(:,1);
	yend=ystart+len*Xhat(:,2);
	hx=quiver(xstart,ystart,xend-xstart,yend-ystart,'Color','blue','ShowArrowHead','on','AutoScale','off');
	%plot Y axis
	xend=xstart+len*Yhat(:,1);
	yend=ystart+len*Yhat(:,2);
	hy=quiver(xstart,ystart,xend-xstart,yend-ystart,'Color','red','ShowArrowHead','on','AutoScale','off');

	legend([hx,hy],'local X direction','local Y direction')
else
	%plot mesh
	A=elements(:,1); B=elements(:,2); C=elements(:,3); D=elements(:,4); E=elements(:,5); F=elements(:,6);
	h1=patch( 'Faces', [A B C],'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
	patch( 'Faces', [D E F],  'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
	patch( 'Faces', [A B E D],'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
	patch( 'Faces', [B E F C],'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
	patch( 'Faces', [C A D F],'Vertices', [x y z],'FaceVertexCData', [1 1 1],'FaceColor','none','EdgeColor','black');
	hold on;

	xstart=x;
	ystart=y;
	zstart=z;
	edgex=max(md.mesh.x(elements),[],2)-min(md.mesh.x(elements),[],2);
	edgez=max(md.mesh.z(elements),[],2)-min(md.mesh.z(elements),[],2);
	len=min(edgex)/1.5;
	lenz=min(edgez)/1.5;
	%plot X axis
	xend=xstart+len*Xhat(:,1);
	yend=ystart+len*Xhat(:,2);
	zend=zstart+len*Xhat(:,3);
	hx=quiver3(xstart,ystart,zstart,xend-xstart,yend-ystart,zend-zstart,'Color','blue','ShowArrowHead','on','AutoScale','off');
	%plot Y axis
	xend=xstart+len*Yhat(:,1);
	yend=ystart+len*Yhat(:,2);
	zend=zstart+len*Yhat(:,3);
	hy=quiver3(xstart,ystart,zstart,xend-xstart,yend-ystart,zend-zstart,'Color','red','ShowArrowHead','on','AutoScale','off');
	%plot Z axis
	xend=xstart+lenz*Zhat(:,1);
	yend=ystart+lenz*Zhat(:,2);
	zend=zstart+lenz*Zhat(:,3);
	hz=quiver3(xstart,ystart,zstart,xend-xstart,yend-ystart,zend-zstart,'Color','green','ShowArrowHead','on','AutoScale','off');

	legend([hx,hy,hz],'local X direction','local Y direction','local Z direction')
end

%apply options
options=addfielddefault(options,'title','Stressbalance referential');
options=addfielddefault(options,'colorbar',0);
applyoptions(md,[],options);
