function [distance] = DistanceToMaskBoundary(x,y,mask);
%DISTANCETOMASKBOUNDARY - Compute distance from any point in a mesh to a mask boundary
%
%   This is a multi-threaded mex file
%   
%   Usage:
%      [distance]=DistanceToMaskBoundary(x,y,mask)
%   
%   x,y,mask: mesh vertices with corresponding mask values.
%   distance: distance from x,y to the mask transition between 0 and 1

% Check usage
if nargin~=3
	help DistanceToMaskBoundary
	error('Wrong usage (see above)');
end

% Call mex module
[distance] = DistanceToMaskBoundary_matlab(x,y,mask);
