function explink(domainoutline,minthreshold,step)
%EXPLINK - allow to link several segments of domain outline together
%
%   Takes a domain outline made of various segments, and links them together in one 
%   domain outline. Use expview to see end result.
%
%   Usage:
%      explink(domainoutline,minthreshold,step)
%
%   See also EXPMASTER, EXPDOC

notdone=1;

while notdone,

	for i=1:1000,
		status=expconcatenate(domainoutline,minthreshold+(i-1)*step);
		if status==0,
			return;
		end
		if status==1,
			break;
		end
		if status==-1,
			continue;
		end
	end
end
