function quality=MeshQuality(md,epsilon,hmin,hmax)
%MESHQUALITY - compute mesh quality
%
%   Usage:
%      MeshQuality(md,epsilon,hmin,hmax);

%Get some variables from the model
index=md.mesh.elements;
x=md.mesh.x;
y=md.mesh.y;

%2d geometric parameter (do not change)
scale=2/9; 

%Compute Hessian
hessian=ComputeHessian(index,x,y,md.inversion.vel_obs,'node');

%Compute metric
metric=ComputeMetric(hessian,scale,epsilon,hmin,hmax,[]);

%Get Areas
areas=GetAreas(index,x,y);

%length edges vectors
e1x=[x(index(:,2))-x(index(:,1))];
e1y=[y(index(:,2))-y(index(:,1))];
e2x=[x(index(:,3))-x(index(:,2))];
e2y=[y(index(:,3))-y(index(:,2))];
e3x=[x(index(:,1))-x(index(:,3))];
e3y=[y(index(:,1))-y(index(:,3))];

%metric of each the 3 nodes for each element
M1=metric(index(:,1),:);
M2=metric(index(:,2),:);
M3=metric(index(:,3),:);

%Get edge length in the metric
L1=1/2*(sqrt(e2x.*(M2(:,1).*e2x+M2(:,2).*e2y)+e2y.*(M2(:,2).*e2x+M2(:,3).*e2y))+sqrt(e1x.*(M1(:,1).*e1x+M1(:,2).*e1y)+e1y.*(M1(:,2).*e1x+M1(:,3).*e1y)));
L2=1/2*(sqrt(e3x.*(M3(:,1).*e3x+M3(:,2).*e3y)+e3y.*(M3(:,2).*e3x+M3(:,3).*e3y))+sqrt(e2x.*(M2(:,1).*e2x+M2(:,2).*e2y)+e2y.*(M2(:,2).*e2x+M2(:,3).*e2y)));
L3=1/2*(sqrt(e1x.*(M1(:,1).*e1x+M1(:,2).*e1y)+e1y.*(M1(:,2).*e1x+M1(:,3).*e1y))+sqrt(e3x.*(M3(:,1).*e3x+M3(:,2).*e3y)+e3y.*(M3(:,2).*e3x+M3(:,3).*e3y)));

%area in the metric
V=1/3*areas.*(sqrt(M1(:,1).*M1(:,3)-M1(:,2).^2)+sqrt(M2(:,1).*M2(:,3)-M2(:,2).^2)+sqrt(M3(:,1).*M3(:,3)-M3(:,2).^2));

%compute quality:
quality=4*sqrt(3)*V./(L1+L2+L3);

%compute error
a=hessian(:,1); b=hessian(:,2); d=hessian(:,3);
a=a(index)*[1;1;1]/3;
b=b(index)*[1;1;1]/3;
d=d(index)*[1;1;1]/3;
lambda1=0.5*((a+d)+sqrt(4*b.^2+(a-d).^2));
lambda2=0.5*((a+d)-sqrt(4*b.^2+(a-d).^2));
lambda1=min(max(abs(lambda1)*scale/epsilon,1/hmax^2),1/hmin^2);
lambda2=min(max(abs(lambda2)*scale/epsilon,1/hmax^2),1/hmin^2);
if length(md.nodeonwater)==md.mesh.numberofvertices;
	pos=find(md.nodeonwater);
	lambda1(pos)=0;
	lambda2(pos)=0;
end
lambda1=lambda1(index)*[1;1;1]/3;
lambda2=lambda2(index)*[1;1;1]/3;

lambdamax=max(lambda1,lambda2);
hmax=max(max(sqrt(e1x.^2+e1y.^2),sqrt(e2x.^2+e2y.^2)),sqrt(e3x.^2+e3y.^2));
epsilon=scale*hmax.^2.*lambdamax;

%display
%X=0:0.1:4; hist(quality,X); xlim([0 3]); title('mesh quality distribution','FontSize',14);
%plotmodel(md,'data',epsilon,'title','Interpolation error','figure',2)
disp(sprintf('\n%s','Mesh Quality'));
disp(sprintf('   %s %g','Average Mesh quality: ',mean(quality)));
disp(sprintf('   %s %g','Worst Element quality:',max(quality)));
disp(sprintf('\n%s','Interpolation Error'));
disp(sprintf('   %s %g %s','Average interpolation error:',mean(epsilon),'m/yr'));
disp(sprintf('   %s %g %s','Maximum interpolation error:',max(epsilon),'m/yr'));
