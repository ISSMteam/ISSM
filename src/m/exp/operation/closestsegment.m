function [profsel indsel]=closestsegment(A,numprofiles,xi,yi)
%CLOSESTSEGMENT - find the closest segment of a profile
%
%   This routine find the segment of the profile A that is the closest
%   to (xi,yi) and return the number of the profile and the number of
%   the first point belonging to this closest segment
%
%   Usage:
%     [profsel indsel]=closestsegment(A,numprofiles,xi,yi) 

	%loop over the middles of each profile, find the closest to (xi,yi)
	profsel=0;
	indsel=0;
	first=1;
	for i=1:numprofiles,
		if length(A(i).x)>1
			middles=[(A(i).x(1:end-1)+A(i).x(2:end))/2 (A(i).y(1:end-1)+A(i).y(2:end))/2];
			distance=(xi-middles(:,1)).^2+(yi-middles(:,2)).^2;
			[newdistance p]=min(distance);
			if (first | (newdistance<olddistance)),
				first=0;
				indsel=p;
				profsel=i;
				olddistance=newdistance;
			end
		end
	end
end
