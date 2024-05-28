function string2=discardnum(string)
%DISCARDNUM -  ??????
%
%   Usage:
%      string2=discardnum(string)

string2=string;

for i=1:length(string),
	if (((string(i)-0) <=57) & ((string(i)-0) >=48)),
		string2=string(1:i-1);
		break;
	end
end
