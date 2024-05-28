function metric=ComputeMetric(hessian,scale,epsilon,hmin,hmax,pos)
%COMPUTEMETRIC - compute metric from an Hessian
%
%   Usage:
%      metric=ComputeMetric(hessian,scale,epsilon,hmin,hmax,pos)
%      pos is contains the positions where the metric is wished to be maximized (water?)
%
%   Example:
%      metric=ComputeMetric(hessian,2/9,10^-1,100,10^5,[])

%first, find the eigen values of each line of H=[hessian(i,1) hessian(i,2); hessian(i,2) hessian(i,3)]
a=hessian(:,1); b=hessian(:,2); d=hessian(:,3);
lambda1=0.5*((a+d)+sqrt(4.*b.^2+(a-d).^2));
lambda2=0.5*((a+d)-sqrt(4.*b.^2+(a-d).^2));
pos1=find(lambda1==0.);
pos2=find(lambda2==0.);
pos3=find(b==0. & lambda1==lambda2);

%Modify the eigen values to control the shape of the elements
lambda1=min(max(abs(lambda1)*scale/epsilon,1./hmax^2),1./hmin^2);
lambda2=min(max(abs(lambda2)*scale/epsilon,1./hmax^2),1./hmin^2);

%compute eigen vectors
norm1=sqrt(8.*b.^2+2.*(d-a).^2+2.*(d-a).*sqrt((a-d).^2+4.*b.^2));
v1x=2.*b./norm1;
v1y=((d-a)+sqrt((a-d).^2+4.*b.^2))./norm1;
norm2=sqrt(8.*b.^2+2.*(d-a).^2-2.*(d-a).*sqrt((a-d).^2+4.*b.^2));
v2x=2.*b./norm2;
v2y=((d-a)-sqrt((a-d).^2+4.*b.^2))./norm2;

v1x(pos3)=1.; v1y(pos3)=0.;
v2x(pos3)=0.; v2y(pos3)=1.;

%Compute new metric (for each node M=V*Lambda*V^-1)
metric=full([(v1x.*v2y-v1y.*v2x).^(-1).*(lambda1.*v2y.*v1x-lambda2.*v1y.*v2x) ...
	(v1x.*v2y-v1y.*v2x).^(-1).*(lambda1.*v1y.*v2y-lambda2.*v1y.*v2y) ...
	(v1x.*v2y-v1y.*v2x).^(-1).*(-lambda1.*v2x.*v1y+lambda2.*v1x.*v2y)]);

%some corrections for 0 eigen values
metric(pos1,:)=repmat([1./hmax^2 0. 1./hmax^2],length(pos1),1);
metric(pos2,:)=repmat([1./hmax^2 0. 1./hmax^2],length(pos2),1);

%take care of water elements
metric(pos,:)=repmat([1./hmax^2 0. 1./hmax^2],length(pos),1);

%take care of NaNs if any (use Matlab eig in a loop)
[pos posj]=find(isnan(metric)); clear posj;
if ~isempty(pos),
	fprintf(' %i %s',length(pos),'NaN found in the metric. Use Matlab routine...');
	for i=1:length(pos)
		H=[hessian(pos(i),1) hessian(pos(i),2)
		hessian(pos(i),2) hessian(pos(i),3)];
		[u,v]=eig(full(H));
		lambda1=v(1,1);
		lambda2=v(2,2);
		v(1,1)=min(max(abs(lambda1)*scale/epsilon,1./hmax^2),1./hmin^2);
		v(2,2)=min(max(abs(lambda2)*scale/epsilon,1./hmax^2),1./hmin^2);

		metricTria=u*v*u^(-1);
		metric(pos(i),:)=[metricTria(1,1) metricTria(1,2) metricTria(2,2)];
	end
end

if any(isnan(metric)),
	error('ComputeMetric error message: NaN in the metric despite our efforts...')
end
