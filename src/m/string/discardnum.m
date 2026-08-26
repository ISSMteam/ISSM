function string2=discardnum(string)
%DISCARDNUM - truncate a string at its first digit
%
%   Return the leading, non-numeric portion of a string, discarding the
%   first digit found and everything that follows it.
%
%   Usage:
%      string2=discardnum(string)

string2=string;

for i=1:length(string)
	if (((string(i)-0) <=57) & ((string(i)-0) >=48))
		string2=string(1:i-1);
		break;
	end
end
