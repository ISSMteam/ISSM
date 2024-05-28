function [xsegments ysegments]=flagedges(elements,x,y,partition)
%FLAGEDGES - return pairs of x,y segments, delimiting partitions.
%
%   Usage:
%      [xsegments ysegments]=flagedges(elements,x,y,partition)

xsegments=[];
ysegments=[];

for i=1:size(elements,1),
	m1=partition(elements(i,1));
	m2=partition(elements(i,2));
	m3=partition(elements(i,3));
	x1=x(elements(i,1));
	x2=x(elements(i,2));
	x3=x(elements(i,3));
	y1=y(elements(i,1));
	y2=y(elements(i,2));
	y3=y(elements(i,3));

	if (m1~=m2) & (m1~=m3) & (m2~=m3),
		xmiddle=(x1+x2+x3)/3;
		ymiddle=(y1+y2+y3)/3;
		xsegments=[xsegments; (x1+x2)/2 xmiddle];
		xsegments=[xsegments; (x1+x3)/2 xmiddle];
		xsegments=[xsegments; (x2+x3)/2 xmiddle];
		ysegments=[ysegments; (y1+y2)/2 ymiddle];
		ysegments=[ysegments; (y1+y3)/2 ymiddle];
		ysegments=[ysegments; (y2+y3)/2 ymiddle];
	end

	if (m1==m2) & (m1~=m3),
		xsegments=[xsegments; (x1+x3)/2 (x2+x3)/2];
		ysegments=[ysegments; (y1+y3)/2 (y2+y3)/2];
	end
	if (m1==m3) & (m2~=m3),
		xsegments=[xsegments; (x1+x2)/2 (x2+x3)/2];
		ysegments=[ysegments; (y1+y2)/2 (y2+y3)/2];
	end

	if (m2==m3) & (m1~=m3),
		xsegments=[xsegments; (x1+x2)/2 (x1+x3)/2];
		ysegments=[ysegments; (y1+y2)/2 (y1+y3)/2];
	end
end
