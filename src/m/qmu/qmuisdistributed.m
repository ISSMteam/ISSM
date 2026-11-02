function found=qmuisdistribted(string)
%QMUISDISTRIBTED - figure out if a string is a descriptor with a numerical postfix, like thickness1, or drag10
%
%   Usage:
%      found=qmuisdistribted(string)

%just take last string element, and see if it is numeric.
last=string(end);

if ((double(last)<=57) & (double(last)>=48))
	found=1;
else
	found=0;
end
