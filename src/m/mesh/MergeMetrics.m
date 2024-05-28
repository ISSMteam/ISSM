function metric=MergeMetrics(metric1,metric2)

M1xx=metric1(:,1); M1xy=metric1(:,2); M1yy=metric1(:,3);
M1lambda1=0.5*((M1xx+M1yy)+sqrt(4*M1xy.^2+(M1xx-M1yy).^2));
M1lambda2=0.5*((M1xx+M1yy)-sqrt(4*M1xy.^2+(M1xx-M1yy).^2));
lambda1=min(M1lambda1,M1lambda2);

M2xx=metric1(:,1); M2xy=metric1(:,2); M2yy=metric1(:,3);
M2lambda1=0.5*((M2xx+M2yy)+sqrt(4*M2xy.^2+(M2xx-M2yy).^2));
M2lambda2=0.5*((M2xx+M2yy)-sqrt(4*M2xy.^2+(M2xx-M2yy).^2));
lambda2=min(M2lambda2,M2lambda2);

metric=metric1;
pos=find(lambda2<lambda1);
metric(pos,:)=metric2(pos,:);
