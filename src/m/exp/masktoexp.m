function masktoexp(x,y,mask,threshold,filename)
%MASKTOEXP - mask to exp file
%
%   Usage:
%      masktoexp(x,y,mask,threshold,filename)
%
%   Example:
%      if A is a matrix of 0 and 1, and we want an exp for
%      the transition:
%      masktoexp(x,y,A,0.5,'contour.exp');
%      To be combined with ExpSimplify

%Create contour for threshold
c=contourc(double(x),double(y),double(mask),[threshold threshold]);
done=0; i=1; j=1;
while (i<length(c))
	num=c(2,i); i=i+1;
	s(j).x=c(1,i:(i+num-1));
	s(j).y=c(2,i:(i+num-1));
	s(j).v=c(1,i);
	i=i+num; j=j+1;
end;

%Create exp structure
A=struct();
if(j-1<1), error('no contour found'); end
for i=1:j-1,
	A(i).x=s(i).x;
	A(i).y=s(i).y;
end;

%write exp
expwrite(A,filename);
