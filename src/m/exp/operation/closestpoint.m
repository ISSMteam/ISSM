function [profsel indsel]=closestpoint(A,numprofiles,xi,yi)
%CLOSESTPOINT - find the closest point of a profile
%
%   This routine find the point of the profile A that is the closest
%   to (xi,yi) and return the number of the profile and the number of
%   the point
%
%   Usage:
%     [profsel indsel]=closestpoint(A,numprofiles,xi,yi) 

	%loop over the points of each profile, find the closest to (xi,yi)
	for i=1:numprofiles,
		distance=(xi-A(i).x).^2+(yi-A(i).y).^2;
		[newdistance p]=min(distance);
		if ((i==1) | (newdistance<olddistance)),
			indsel=p;
			profsel=i;
			olddistance=newdistance;
		end
	end
end
