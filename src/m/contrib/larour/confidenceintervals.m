function level=confidenceintervals(x,h,threshold)
%CONFIDENCEINTERVALS: build confidence interval from histogram.
   
	nods=size(h,1); nbins=size(h,2);

	%make sure h is one column smaller than x
	if size(x,2) ~= (size(h,2)+1),
		error('x should be one column larger than the histogram');
	end

	%First build cdf: 
	c=[zeros(nods,1) cumsum(h,2)];

	%In percentage: 
	threshold=threshold/100;

	flags=zeros(nods,nbins+1);
	pos=find(c<threshold);
	flags(pos)=1;
	diffc=diff(flags,1,2);

	[ind,j]=find(diffc==-1);
	indices=zeros(nods,1); 
	indices(ind)=j;

	%threshold was passed between j and j+1 for each node. try 
	%to pin down the fraction better: 
	idx = sub2ind(size(c), [1:size(c,1)]', indices);
	idx2 = sub2ind(size(c), [1:size(c,1)]', indices+1);
	val1=c(idx);
	val2=c(idx2);
	x1=x(idx);
	x2=x(idx2);

	fraction=(threshold-val1)./(val2-val1);
	level=x1+fraction.*(x2-x1);
